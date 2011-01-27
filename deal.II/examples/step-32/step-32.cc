//TODO: - adjust stopping criteria for solvers
//      - better refinement at the start?
//      - check solver stability
//      - Q2 Mapping useful?


/* $Id$ */
/* Author: Martin Kronbichler, Uppsala University,
           Wolfgang Bangerth, Texas A&M University 2008, 2009 */
/*                                                                */
/*    Copyright (C) 2008, 2009, 2010, 2011 by the deal.II authors */
/*                                                                */
/*    This file is subject to QPL and may not be  distributed     */
/*    without copyright and license information. Please refer     */
/*    to the file deal.II/doc/license.html for the  text  and     */
/*    further information on this license.                        */

				 // @sect3{Include files}

				 // We include the functionality
				 // of these well-known deal.II
				 // library files and some C++
				 // header files.
#include <base/quadrature_lib.h>
#include <base/logstream.h>
#include <base/function.h>
#include <base/utilities.h>
#include <base/conditional_ostream.h>
#include <base/work_stream.h>
#include <base/timer.h>

#include <lac/full_matrix.h>
#include <lac/solver_bicgstab.h>
#include <lac/solver_cg.h>
#include <lac/solver_gmres.h>
#include <lac/constraint_matrix.h>
#include <lac/block_sparsity_pattern.h>
#include <lac/trilinos_block_vector.h>
#include <lac/trilinos_sparse_matrix.h>
#include <lac/trilinos_block_sparse_matrix.h>
#include <lac/trilinos_precondition.h>
#include <lac/trilinos_solver.h>

#include <grid/tria.h>
#include <grid/grid_generator.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/filtered_iterator.h>
#include <grid/tria_boundary_lib.h>
#include <grid/grid_tools.h>
#include <grid/grid_refinement.h>

#include <dofs/dof_handler.h>
#include <dofs/dof_renumbering.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>

#include <fe/fe_q.h>
#include <fe/fe_dgq.h>
#include <fe/fe_system.h>
#include <fe/fe_values.h>
#include <fe/mapping_q1.h>

#include <numerics/vectors.h>
#include <numerics/matrices.h>
#include <numerics/data_out.h>
#include <numerics/error_estimator.h>
#include <numerics/solution_transfer.h>

#include <fstream>
#include <iostream>
#include <limits>

				 // This is the only include file that is new:
				 // We use an IndexSet to describe the
				 // %parallel partitioning of vectors and
				 // matrices.
#include <base/index_set.h>

#include <distributed/tria.h>
#include <distributed/solution_transfer.h>
#include <distributed/grid_refinement.h>

#include <iostream>
#include <sstream>
#include <string>



using namespace dealii;

				 // @sect3{Equation data}

				 // In the following namespace, we define the
				 // various pieces of equation data. All of
				 // these are exhaustively discussed in the
				 // description of the testcase in the
				 // introduction:
namespace EquationData
{
  const double eta                   = 1e21;    /* Pa s       */
  const double kappa                 = 1e-6;
  const double reference_density     = 3300;    /* kg / m^3   */
  const double reference_temperature = 293;     /* K          */
  const double expansion_coefficient = 2e-5;    /* 1/K        */
  const double specific_heat         = 1250;    /* J / K / kg */  //??
  const double radiogenic_heating    = 7.4e-12; /* W / kg     */  //??

  const double R0      = 6371000.-2890000.;     /* m          */
  const double R1      = 6371000.-  35000.;     /* m          */

  const double T0      = 4000+273;              /* K          */
  const double T1      =  700+273;              /* K          */

  const double year_in_seconds  = 60*60*24*365.2425;
  const double end_time         = 1e8 * year_in_seconds;

  const double pressure_scaling = eta / (R1-R0);


  double density (const double temperature)
  {
    return (reference_density *
	    (1 - expansion_coefficient * (temperature -
					  reference_temperature)));
  }


  template <int dim>
  Tensor<1,dim> gravity_vector (const Point<dim> &p)
  {
// interpolate the following values with a physically realistic model:
//    const double g0      = 10.7;                  /* m / s^2    */
//    const double g1      = 9.81;                  /* m / s^2    */
  
    const double r = p.norm();
    return -(1.245e-6 * r + 7.714e13/r/r) * p / p.norm();
  }


  template <int dim>
  class TemperatureInitialValues : public Function<dim>
  {
    public:
      TemperatureInitialValues () : Function<dim>(1) {}

      virtual double value (const Point<dim>   &p,
			    const unsigned int  component = 0) const;

      virtual void vector_value (const Point<dim> &p,
				 Vector<double>   &value) const;
  };



  template <int dim>
  double
  TemperatureInitialValues<dim>::value (const Point<dim>  &p,
					const unsigned int) const
  {
    const double r = p.norm();
    const double h = R1-R0;

				     // s = fraction of the way from
				     // the inner to the outer
				     // boundary; 0<=s<=1
    const double s = (r-R0)/h;

/* now compute an angular variation of the linear temperature field by
   stretching the variable s appropriately. note that the following
   formula leaves the end points s=0 and s=1 fixed, but stretches the
   region in between depending on the angle phi=atan2(x,y).

   For a plot, see
   http://www.wolframalpha.com/input/?i=plot+%28%282*sqrt%28x^2%2By^2%29-1%29%2B0.2*%282*sqrt%28x^2%2By^2%29-1%29*%281-%282*sqrt%28x^2%2By^2%29-1%29%29*sin%286*atan2%28x%2Cy%29%29%29%2C+x%3D-1+to+1%2C+y%3D-1+to+1
*/
    const double phi   = std::atan2(p(0),p(1));
    const double s_mod = (s
			  +
			  0.2 * s * (1-s) * std::sin(6*phi));

    return T0*(1.0-s_mod) + T1*s_mod;
  }


  template <int dim>
  void
  TemperatureInitialValues<dim>::vector_value (const Point<dim> &p,
					       Vector<double>   &values) const
  {
    for (unsigned int c=0; c<this->n_components; ++c)
      values(c) = TemperatureInitialValues<dim>::value (p, c);
  }
}



				 // @sect3{Linear solvers and preconditioners}

				 // In comparison to step-31, we did one
				 // change in the linear algebra of the
				 // problem: We exchange the
				 // <code>InverseMatrix</code> that
				 // previously held the approximation of the
				 // Schur complement by a preconditioner
				 // only (we will choose ILU in the
				 // application code below), as discussed in
				 // the introduction. This trick we already
				 // did for the velocity block - the idea of
				 // this is that the solver iterations on
				 // the block system will eventually also
				 // make the approximation for the Schur
				 // complement good. If the preconditioner
				 // we're using is good enough, there will
				 // be no increase in the outer iteration
				 // count compared to using converged solves
				 // for the inverse matrices of velocity and
				 // Schur complement. All we need to do for
				 // implementing that change is to give the
				 // respective variable in the
				 // BlockSchurPreconditioner class another
				 // name.
namespace LinearSolvers
{
  template <class PreconditionerA, class PreconditionerMp>
  class RightPrecond : public Subscriptor
  {
    public:
      RightPrecond (
	const TrilinosWrappers::BlockSparseMatrix  &S,
	const TrilinosWrappers::BlockSparseMatrix  &Spre,
	const PreconditionerMp                     &Mppreconditioner,
	const PreconditionerA                      &Apreconditioner)
		  :
		  stokes_matrix     (&S),
		  stokes_preconditioner_matrix     (&Spre),
		  mp_preconditioner (Mppreconditioner),
		  a_preconditioner  (Apreconditioner)
	{}

     void solve_S(TrilinosWrappers::MPI::Vector &dst,
		  const TrilinosWrappers::MPI::Vector &src) const
	{
	  SolverControl cn(5000, 1e-5);//src.l2_norm()*1e-5);

	  TrilinosWrappers::SolverCG solver(cn);

	  solver.solve(stokes_preconditioner_matrix->block(1,1),
		       dst, src,
		       mp_preconditioner);

	  dst*=-1.0;
	}

      void solve_A(TrilinosWrappers::MPI::Vector &dst,
		  const TrilinosWrappers::MPI::Vector &src) const
	{
	  SolverControl cn(5000, src.l2_norm()*1e-4);
	  TrilinosWrappers::SolverCG solver(cn);
	  solver.solve(stokes_matrix->block(0,0), dst, src, a_preconditioner);
	}

      void vmult (TrilinosWrappers::MPI::BlockVector       &dst,
		  const TrilinosWrappers::MPI::BlockVector &src) const
	{
	  TrilinosWrappers::MPI::Vector utmp(src.block(0));

	  solve_S(dst.block(1), src.block(1));

	  stokes_matrix->block(0,1).vmult(utmp, dst.block(1)); //B^T
	  utmp*=-1.0;
	  utmp.add(src.block(0));

	  solve_A(dst.block(0), utmp);
	}

    private:
      const SmartPointer<const TrilinosWrappers::BlockSparseMatrix> stokes_matrix;
      const SmartPointer<const TrilinosWrappers::BlockSparseMatrix> stokes_preconditioner_matrix;
      const PreconditionerMp &mp_preconditioner;
      const PreconditionerA  &a_preconditioner;
  };
}



				 // @sect3{Definition of assembly data structures}
				 //
				 // As described in the introduction, we will
				 // use the WorkStream mechanism discussed in
				 // the @ref threads module to parallelize
				 // operations among the processors of a
				 // single machine. The WorkStream class
				 // requires that data is passed around in two
				 // kinds of data structures, one for scratch
				 // data and one to pass data from the
				 // assembly function to the function that
				 // copies local contributions into global
				 // objects.
				 //
				 // The following namespace (and the two
				 // sub-namespaces) contains a collection of
				 // data structures that serve this purpose,
				 // one pair for each of the four operations
				 // discussed in the introduction that we will
				 // want to parallelize. Each
				 // assembly routine gets two sets of data: a
				 // Scratch array that collects all the
				 // classes and arrays that are used for the
				 // calculation of the cell contribution, and
				 // a CopyData array that keeps local matrices
				 // and vectors which will be written into the
				 // global matrix. Whereas CopyData is a
				 // container for the final data that is
				 // written into the global matrices and
				 // vector (and, thus, absolutely necessary),
				 // the Scratch arrays are merely there for
				 // performance reasons &mdash; it would be
				 // much more expensive to set up a FEValues
				 // object on each cell, than creating it only
				 // once and updating some derivative data.
				 //
				 // Using the program in step-31, we have
				 // four assembly routines. One for the
				 // preconditioner matrix of the Stokes
				 // system, one for the Stokes matrix and
				 // right hand side, one for the
				 // temperature matrices and one for the
				 // right hand side of the temperature
				 // equation. We organize the scratch
				 // arrays and a CopyData arrays for each
				 // of those four assembly components
				 // using a <code>struct</code>
				 // environment.
				 //
				 // Regarding the Scratch array, each
				 // struct is equipped with a constructor
				 // that create an FEValues object for a
				 // @ref FiniteElement "finite element", a
				 // @ref Quadrature "quadrature formula"
				 // and some
				 // @ref UpdateFlags "update flags".
				 // Moreover, we manually
				 // implement a copy constructor (since
				 // the FEValues class is not copyable by
				 // itself), and provide some additional
				 // vector fields that are used to improve
				 // performance of assembly.
namespace Assembly
{
  namespace Scratch
  {
    template <int dim>
    struct StokesPreconditioner
    {
	StokesPreconditioner (const FiniteElement<dim> &stokes_fe,
			      const Quadrature<dim>    &stokes_quadrature,
			      const UpdateFlags         update_flags);
	StokesPreconditioner (const StokesPreconditioner &data);

	FEValues<dim>               stokes_fe_values;

	std::vector<Tensor<2,dim> > grad_phi_u;
	std::vector<double>         phi_p;
    };

    template <int dim>
    StokesPreconditioner<dim>::
    StokesPreconditioner (const FiniteElement<dim> &stokes_fe,
			  const Quadrature<dim>    &stokes_quadrature,
			  const UpdateFlags         update_flags)
		    :
		    stokes_fe_values (stokes_fe, stokes_quadrature,
				      update_flags),
		    grad_phi_u (stokes_fe.dofs_per_cell),
		    phi_p (stokes_fe.dofs_per_cell)
    {}



    template <int dim>
    StokesPreconditioner<dim>::
    StokesPreconditioner (const StokesPreconditioner &scratch)
		    :
		    stokes_fe_values (scratch.stokes_fe_values.get_fe(),
				      scratch.stokes_fe_values.get_quadrature(),
				      scratch.stokes_fe_values.get_update_flags()),
		    grad_phi_u (scratch.grad_phi_u),
		    phi_p (scratch.phi_p)
    {}



				     // Observe that we derive the
				     // StokesSystem scratch array from the
				     // StokesPreconditioner array. We do this
				     // because all the objects that are
				     // necessary for the assembly of the
				     // preconditioner are also needed for the
				     // actual matrix system and right hand
				     // side, plus some extra data. This makes
				     // the program more compact. Note also
				     // that the assembly of the Stokes system
				     // and the temperature right hand side
				     // further down requires data from
				     // temperature and velocity,
				     // respectively, so we actually need two
				     // FEValues objects for those two cases.
    template <int dim>
    struct StokesSystem : public StokesPreconditioner<dim>
    {
	StokesSystem (const FiniteElement<dim> &stokes_fe,
		      const Quadrature<dim>    &stokes_quadrature,
		      const UpdateFlags         stokes_update_flags,
		      const FiniteElement<dim> &temperature_fe,
		      const UpdateFlags         temperature_update_flags);

	StokesSystem (const StokesSystem<dim> &data);

	FEValues<dim>  temperature_fe_values;

	std::vector<Tensor<1,dim> >          phi_u;
	std::vector<SymmetricTensor<2,dim> > grads_phi_u;
	std::vector<double>                  div_phi_u;

	std::vector<double>                  old_temperature_values;
    };


    template <int dim>
    StokesSystem<dim>::
    StokesSystem (const FiniteElement<dim> &stokes_fe,
		  const Quadrature<dim>    &stokes_quadrature,
		  const UpdateFlags         stokes_update_flags,
		  const FiniteElement<dim> &temperature_fe,
		  const UpdateFlags         temperature_update_flags)
		    :
		    StokesPreconditioner<dim> (stokes_fe, stokes_quadrature,
					       stokes_update_flags),
		    temperature_fe_values (temperature_fe, stokes_quadrature,
					   temperature_update_flags),
		    phi_u (stokes_fe.dofs_per_cell),
		    grads_phi_u (stokes_fe.dofs_per_cell),
		    div_phi_u (stokes_fe.dofs_per_cell),
		    old_temperature_values (stokes_quadrature.size())
    {}


    template <int dim>
    StokesSystem<dim>::
    StokesSystem (const StokesSystem<dim> &scratch)
		    :
		    StokesPreconditioner<dim> (scratch),
		    temperature_fe_values (scratch.temperature_fe_values.get_fe(),
					   scratch.temperature_fe_values.get_quadrature(),
					   scratch.temperature_fe_values.get_update_flags()),
		    phi_u (scratch.phi_u),
		    grads_phi_u (scratch.grads_phi_u),
		    div_phi_u (scratch.div_phi_u),
		    old_temperature_values (scratch.old_temperature_values)
    {}



    template <int dim>
    struct TemperatureMatrix
    {
	TemperatureMatrix (const FiniteElement<dim> &temperature_fe,
			   const Quadrature<dim>    &temperature_quadrature);
	TemperatureMatrix (const TemperatureMatrix &data);

	FEValues<dim>               temperature_fe_values;

	std::vector<double>         phi_T;
	std::vector<Tensor<1,dim> > grad_phi_T;
    };

    template <int dim>
    TemperatureMatrix<dim>::
    TemperatureMatrix (const FiniteElement<dim> &temperature_fe,
		       const Quadrature<dim>    &temperature_quadrature)
		    :
		    temperature_fe_values (temperature_fe, temperature_quadrature,
					   update_values    | update_gradients |
					   update_JxW_values),
		    phi_T (temperature_fe.dofs_per_cell),
		    grad_phi_T (temperature_fe.dofs_per_cell)
    {}


    template <int dim>
    TemperatureMatrix<dim>::
    TemperatureMatrix (const TemperatureMatrix &scratch)
		    :
		    temperature_fe_values (scratch.temperature_fe_values.get_fe(),
					   scratch.temperature_fe_values.get_quadrature(),
					   scratch.temperature_fe_values.get_update_flags()),
		    phi_T (scratch.phi_T),
		    grad_phi_T (scratch.grad_phi_T)
    {}


    template <int dim>
    struct TemperatureRHS
    {
	TemperatureRHS (const FiniteElement<dim> &temperature_fe,
			const FiniteElement<dim> &stokes_fe,
			const Quadrature<dim>    &quadrature);
	TemperatureRHS (const TemperatureRHS &data);

	FEValues<dim>               temperature_fe_values;
	FEValues<dim>               stokes_fe_values;

	std::vector<double>         phi_T;
	std::vector<Tensor<1,dim> > grad_phi_T;

	std::vector<Tensor<1,dim> > old_velocity_values;
	std::vector<Tensor<1,dim> > old_old_velocity_values;

	std::vector<SymmetricTensor<2,dim> > old_strain_rates;
	std::vector<SymmetricTensor<2,dim> > old_old_strain_rates;
	
	std::vector<double>         old_temperature_values;
	std::vector<double>         old_old_temperature_values;
	std::vector<Tensor<1,dim> > old_temperature_grads;
	std::vector<Tensor<1,dim> > old_old_temperature_grads;
	std::vector<double>         old_temperature_laplacians;
	std::vector<double>         old_old_temperature_laplacians;
    };

    template <int dim>
    TemperatureRHS<dim>::
    TemperatureRHS (const FiniteElement<dim> &temperature_fe,
		    const FiniteElement<dim> &stokes_fe,
		    const Quadrature<dim>    &quadrature)
		    :
		    temperature_fe_values (temperature_fe, quadrature,
					   update_values    |
					   update_gradients |
					   update_hessians  |
					   update_quadrature_points |
					   update_JxW_values),
		    stokes_fe_values (stokes_fe, quadrature,
				      update_values | update_gradients),
		    phi_T (temperature_fe.dofs_per_cell),
		    grad_phi_T (temperature_fe.dofs_per_cell),

		    old_velocity_values (quadrature.size()),
		    old_old_velocity_values (quadrature.size()),
		    old_strain_rates (quadrature.size()),
		    old_old_strain_rates (quadrature.size()),
		    
		    old_temperature_values (quadrature.size()),
		    old_old_temperature_values(quadrature.size()),
		    old_temperature_grads(quadrature.size()),
		    old_old_temperature_grads(quadrature.size()),
		    old_temperature_laplacians(quadrature.size()),
		    old_old_temperature_laplacians(quadrature.size())
    {}


    template <int dim>
    TemperatureRHS<dim>::
    TemperatureRHS (const TemperatureRHS &scratch)
		    :
		    temperature_fe_values (scratch.temperature_fe_values.get_fe(),
					   scratch.temperature_fe_values.get_quadrature(),
					   scratch.temperature_fe_values.get_update_flags()),
		    stokes_fe_values (scratch.stokes_fe_values.get_fe(),
				      scratch.stokes_fe_values.get_quadrature(),
				      scratch.stokes_fe_values.get_update_flags()),
		    phi_T (scratch.phi_T),
		    grad_phi_T (scratch.grad_phi_T),

		    old_velocity_values (scratch.old_velocity_values),
		    old_old_velocity_values (scratch.old_old_velocity_values),
		    old_strain_rates (scratch.old_strain_rates),
		    old_old_strain_rates (scratch.old_old_strain_rates),

		    old_temperature_values (scratch.old_temperature_values),
		    old_old_temperature_values (scratch.old_old_temperature_values),
		    old_temperature_grads (scratch.old_temperature_grads),
		    old_old_temperature_grads (scratch.old_old_temperature_grads),
		    old_temperature_laplacians (scratch.old_temperature_laplacians),
		    old_old_temperature_laplacians (scratch.old_old_temperature_laplacians)
    {}
  }


				   // The CopyData arrays are similar to the
				   // Scratch arrays. They provide a
				   // constructor, a copy operation, and
				   // some arrays for local matrix, local
				   // vectors and the relation between local
				   // and global degrees of freedom (a.k.a.
				   // <code>local_dof_indices</code>).
  namespace CopyData
  {
    template <int dim>
    struct StokesPreconditioner
    {
	StokesPreconditioner (const FiniteElement<dim> &stokes_fe);
	StokesPreconditioner (const StokesPreconditioner &data);

	FullMatrix<double>          local_matrix;
	std::vector<unsigned int>   local_dof_indices;
    };

    template <int dim>
    StokesPreconditioner<dim>::
    StokesPreconditioner (const FiniteElement<dim> &stokes_fe)
		    :
		    local_matrix (stokes_fe.dofs_per_cell,
				  stokes_fe.dofs_per_cell),
		    local_dof_indices (stokes_fe.dofs_per_cell)
    {}



    template <int dim>
    StokesPreconditioner<dim>::
    StokesPreconditioner (const StokesPreconditioner &data)
		    :
		    local_matrix (data.local_matrix),
		    local_dof_indices (data.local_dof_indices)
    {}



    template <int dim>
    struct StokesSystem : public StokesPreconditioner<dim>
    {
	StokesSystem (const FiniteElement<dim> &stokes_fe);
	StokesSystem (const StokesSystem<dim> &data);

	Vector<double> local_rhs;
    };


    template <int dim>
    StokesSystem<dim>::
    StokesSystem (const FiniteElement<dim> &stokes_fe)
		    :
		    StokesPreconditioner<dim> (stokes_fe),
		    local_rhs (stokes_fe.dofs_per_cell)
    {}


    template <int dim>
    StokesSystem<dim>::
    StokesSystem (const StokesSystem<dim> &data)
		    :
		    StokesPreconditioner<dim> (data),
		    local_rhs (data.local_rhs)
    {}



    template <int dim>
    struct TemperatureMatrix
    {
	TemperatureMatrix (const FiniteElement<dim> &temperature_fe);
	TemperatureMatrix (const TemperatureMatrix &data);

	FullMatrix<double>          local_mass_matrix;
	FullMatrix<double>          local_stiffness_matrix;
	std::vector<unsigned int>   local_dof_indices;
    };

    template <int dim>
    TemperatureMatrix<dim>::
    TemperatureMatrix (const FiniteElement<dim> &temperature_fe)
		    :
		    local_mass_matrix (temperature_fe.dofs_per_cell,
				       temperature_fe.dofs_per_cell),
		    local_stiffness_matrix (temperature_fe.dofs_per_cell,
					    temperature_fe.dofs_per_cell),
		    local_dof_indices (temperature_fe.dofs_per_cell)
    {}


    template <int dim>
    TemperatureMatrix<dim>::
    TemperatureMatrix (const TemperatureMatrix &data)
		    :
		    local_mass_matrix (data.local_mass_matrix),
		    local_stiffness_matrix (data.local_stiffness_matrix),
		    local_dof_indices (data.local_dof_indices)
    {}


    template <int dim>
    struct TemperatureRHS
    {
	TemperatureRHS (const FiniteElement<dim> &temperature_fe);
	TemperatureRHS (const TemperatureRHS &data);

	Vector<double>              local_rhs;
	std::vector<unsigned int>   local_dof_indices;
        FullMatrix<double>          matrix_for_bc;
    };

    template <int dim>
    TemperatureRHS<dim>::
    TemperatureRHS (const FiniteElement<dim> &temperature_fe)
		    :
		    local_rhs (temperature_fe.dofs_per_cell),
		    local_dof_indices (temperature_fe.dofs_per_cell),
		    matrix_for_bc (temperature_fe.dofs_per_cell,
				   temperature_fe.dofs_per_cell)
    {}


    template <int dim>
    TemperatureRHS<dim>::
    TemperatureRHS (const TemperatureRHS &data)
		    :
		    local_rhs (data.local_rhs),
		    local_dof_indices (data.local_dof_indices),
		    matrix_for_bc (data.matrix_for_bc)
    {}
  }
}



				 // @sect3{The <code>BoussinesqFlowProblem</code> class template}
				 //
				 // This is the declaration of the main
				 // class. It is very similar to
				 // step-31. Following the @ref
				 // MTWorkStream "task-based parallelization"
				 // paradigm, we split all the
				 // assembly routines into two parts: a
				 // first part that can do all the
				 // calculations on a certain cell without
				 // taking care of other threads, and a
				 // second part (which is writing the
				 // local data into the global matrices
				 // and vectors) which can be entered by
				 // only one thread at a time. In order to
				 // implement that, we provide functions
				 // for each of those two steps for all
				 // the four assembly routines that we use
				 // in this program.
				 //
				 // The <code>pcout</code> (for <i>%parallel
				 // <code>std::cout</code></i>) object is used
				 // to simplify writing output: each MPI
				 // process can use this to generate output as
				 // usual, but since each of these processes
				 // will produce the same output it will just
				 // be replicated many times over; with the
				 // ConditionalOStream class, only the output
				 // generated by one task will actually be
				 // printed to screen, whereas the output by
				 // all the other threads will simply be
				 // forgotten.
				 //
				 // In a bit of naming confusion, you will
				 // notice below that some of the variables
				 // from namespace TrilinosWrappers are
				 // taken from namespace
				 // TrilinosWrappers::MPI (such as the right
				 // hand side vectors) whereas others are
				 // not (such as the various matrices). For
				 // the matrices, we happen to use the same
				 // class names for %parallel and sequential
				 // data structures, i.e. all matrices will
				 // actually be considered %parallel
				 // below. On the other hand, for vectors,
				 // only those from namespace
				 // TrilinosWrappers::MPI are actually
				 // distributed. In particular, we will
				 // frequently have to query velocities and
				 // temperatures at arbitrary quadrature
				 // points; consequently, rather than
				 // "localizing" a vector whenever we need a
				 // localized vector, we solve linear
				 // systems in %parallel but then immediately
				 // localize the solution for further
				 // processing. The various
				 // <code>*_solution</code> vectors are
				 // therefore filled immediately after
				 // solving their respective linear system
				 // in %parallel.
				 //
				 // The only other new data member is
				 // <code>computing_timer</code>. Its class
				 // type, TimerOutput, can be used to
				 // conveniently account for compute time
				 // spent in certain "sections" of the code
 				 // that are repeatedly entered. For
 				 // example, we will enter (and leave)
 				 // sections for Stokes matrix assembly and
 				 // would like to accumulate the run time
 				 // spent in this section over all time
 				 // steps. At the end of the program, the
 				 // destructor of the TimerOutput class will
 				 // automatically produce a nice summary of
 				 // the times spent in all the sections. For
 				 // this output, one can choose whether wall
 				 // clock or CPU times are to be printed, as
 				 // well as whether we want to produce
				 // output every time we leave a section --
 				 // which would be quite a lot of additional
 				 // output -- or just in the end of the
 				 // program (this choice is made in the
 				 // from this variable in the results
 				 // section of this tutorial program.
template <int dim>
class BoussinesqFlowProblem
{
  public:
    BoussinesqFlowProblem ();
    void run (unsigned int ref);

  private:
    void setup_dofs ();
    void assemble_stokes_preconditioner ();
    void build_stokes_preconditioner ();
    void assemble_stokes_system ();
    void assemble_temperature_matrix ();
    void assemble_temperature_system (const double maximal_velocity);
    void project_temperature_field ();
    double get_maximal_velocity () const;
    std::pair<double,double> get_extrapolated_temperature_range () const;
    void solve ();
    void output_results ();
    void refine_mesh (const unsigned int max_grid_level);

    double
    compute_viscosity(const std::vector<double>          &old_temperature,
		      const std::vector<double>          &old_old_temperature,
		      const std::vector<Tensor<1,dim> >  &old_temperature_grads,
		      const std::vector<Tensor<1,dim> >  &old_old_temperature_grads,
		      const std::vector<double>          &old_temperature_laplacians,
		      const std::vector<double>          &old_old_temperature_laplacians,
		      const std::vector<Tensor<1,dim> >  &old_velocity_values,
		      const std::vector<Tensor<1,dim> >  &old_old_velocity_values,
		      const std::vector<SymmetricTensor<2,dim> >  &old_strain_rates,
		      const std::vector<SymmetricTensor<2,dim> >  &old_old_strain_rates,
		      const double                        global_u_infty,
		      const double                        global_T_variation,
		      const double                        cell_diameter) const;


    ConditionalOStream                  pcout;

    parallel::distributed::Triangulation<dim> triangulation;
    double                              global_Omega_diameter;

    const unsigned int                  stokes_degree;
    FESystem<dim>                       stokes_fe;
    DoFHandler<dim> stokes_dof_handler;
    ConstraintMatrix                    stokes_constraints;

    TrilinosWrappers::BlockSparseMatrix stokes_matrix;
    TrilinosWrappers::BlockSparseMatrix stokes_preconditioner_matrix;

    TrilinosWrappers::MPI::BlockVector  stokes_solution;
    TrilinosWrappers::MPI::BlockVector  old_stokes_solution;
    TrilinosWrappers::MPI::BlockVector  stokes_rhs;


    const unsigned int                  temperature_degree;
    FE_Q<dim>                           temperature_fe;
    DoFHandler<dim>                     temperature_dof_handler;
    ConstraintMatrix                    temperature_constraints;

    TrilinosWrappers::SparseMatrix      temperature_mass_matrix;
    TrilinosWrappers::SparseMatrix      temperature_stiffness_matrix;
    TrilinosWrappers::SparseMatrix      temperature_matrix;

    TrilinosWrappers::MPI::Vector       temperature_solution;
    TrilinosWrappers::MPI::Vector       old_temperature_solution;
    TrilinosWrappers::MPI::Vector       old_old_temperature_solution;
    TrilinosWrappers::MPI::Vector       temperature_rhs;


    double time_step;
    double old_time_step;
    unsigned int timestep_number;

    std_cxx1x::shared_ptr<TrilinosWrappers::PreconditionAMG> Amg_preconditioner;
    std_cxx1x::shared_ptr<TrilinosWrappers::PreconditionILU> Mp_preconditioner;
    std_cxx1x::shared_ptr<TrilinosWrappers::PreconditionIC>  T_preconditioner;

    bool rebuild_stokes_matrix;
    bool rebuild_stokes_preconditioner;
    bool rebuild_temperature_matrices;
    bool rebuild_temperature_preconditioner;

    TimerOutput computing_timer;

    void setup_stokes_matrix (const std::vector<IndexSet> &stokes_partitioning);
    void setup_stokes_preconditioner (const std::vector<IndexSet> &stokes_partitioning);
    void setup_temperature_matrices (const IndexSet &temperature_partitioning);

    void
    local_assemble_stokes_preconditioner (const typename DoFHandler<dim>::active_cell_iterator &cell,
					  Assembly::Scratch::StokesPreconditioner<dim> &scratch,
					  Assembly::CopyData::StokesPreconditioner<dim> &data);

    void
    copy_local_to_global_stokes_preconditioner (const Assembly::CopyData::StokesPreconditioner<dim> &data);


    void
    local_assemble_stokes_system (const typename DoFHandler<dim>::active_cell_iterator &cell,
				  Assembly::Scratch::StokesSystem<dim>  &scratch,
				  Assembly::CopyData::StokesSystem<dim> &data);

    void
    copy_local_to_global_stokes_system (const Assembly::CopyData::StokesSystem<dim> &data);


    void
    local_assemble_temperature_matrix (const typename DoFHandler<dim>::active_cell_iterator &cell,
				       Assembly::Scratch::TemperatureMatrix<dim>  &scratch,
				       Assembly::CopyData::TemperatureMatrix<dim> &data);

    void
    copy_local_to_global_temperature_matrix (const Assembly::CopyData::TemperatureMatrix<dim> &data);



    void
    local_assemble_temperature_rhs (const std::pair<double,double> global_T_range,
				    const double                   global_max_velocity,
				    const typename DoFHandler<dim>::active_cell_iterator &cell,
				    Assembly::Scratch::TemperatureRHS<dim> &scratch,
				    Assembly::CopyData::TemperatureRHS<dim> &data);

    void
    copy_local_to_global_temperature_rhs (const Assembly::CopyData::TemperatureRHS<dim> &data);

    class Postprocessor;
};


				 // @sect3{BoussinesqFlowProblem class implementation}

				 // @sect4{BoussinesqFlowProblem::BoussinesqFlowProblem}
				 //
				 // The constructor of the problem is very
				 // similar to the constructor in
				 // step-31. What is different is the
				 // %parallel communication: Trilinos uses a
				 // message passing interface (MPI) for data
				 // distribution. When entering the
				 // BoussinesqFlowProblem class, we have to
				 // decide how the parallization is to be
				 // done. We choose a rather simple strategy
				 // and let all processors that are running
				 // the program work together, specified by
				 // the communicator
				 // <code>comm_world()</code>. Next, we
				 // create some modified output stream as we
				 // already did in step-18. In MPI, all the
				 // processors run the same program
				 // individually (they simply operate on
				 // different chunks of data and exchange
				 // some part of that data from time to
				 // time). Next, we need to initialize the
 				 // <code>pcout</code> object in order to
 				 // print the user information only on one
 				 // processor. The implementation of this
 				 // idea is to check the process number when
 				 // <code>pcout</code> gets a true argument,
 				 // and it uses the <code>std::cout</code>
 				 // stream for output. If we are one
 				 // processor five, for instance, then we
 				 // will give a <code>false</code> argument
 				 // to <code>pcout</code>, which means that
 				 // the output of that processor will not be
 				 // printed anywhere.
 				 //
 				 // Finally, we enter the preferred options
 				 // for the TimerOutput object to its
 				 // constructor. We restrict the output to
 				 // the <code>pcout</code> stream (processor
 				 // 0), and then we specify that we want to
 				 // get a summary table in the end of the
 				 // program which shows us wallclock times
 				 // (as opposed to CPU times).
template <int dim>
BoussinesqFlowProblem<dim>::BoussinesqFlowProblem ()
                :
		pcout (std::cout,
		       (Utilities::System::
			get_this_mpi_process(MPI_COMM_WORLD)
			== 0)),

		triangulation (MPI_COMM_WORLD,
			       typename Triangulation<dim>::MeshSmoothing
			       (Triangulation<dim>::smoothing_on_refinement | Triangulation<dim>::smoothing_on_coarsening
			       )
		),

                stokes_degree (1),
                stokes_fe (FE_Q<dim>(stokes_degree+1), dim,
			   FE_Q<dim>(stokes_degree), 1),
		stokes_dof_handler (triangulation),

		temperature_degree (2),
		temperature_fe (temperature_degree),
                temperature_dof_handler (triangulation),

                time_step (0),
		old_time_step (0),
		timestep_number (0),
		rebuild_stokes_matrix (true),
		rebuild_stokes_preconditioner (true),
		rebuild_temperature_matrices (true),
		rebuild_temperature_preconditioner (true),

		computing_timer (pcout, TimerOutput::summary,
				 TimerOutput::wall_times)
{}



				 // @sect4{The BoussinesqFlowProblem helper functions}
				 //
				 // Except two small details, this
				 // function is the very same as in
				 // step-31. The first detail is
				 // actually common to all functions
				 // that implement loop over all cells
				 // in the triangulation: When
				 // operating in %parallel, each
				 // processor only works on a chunk of
				 // cells. This chunk of cells is
				 // identified via a so-called
				 // <code>subdomain_id</code>, as we
				 // also did in step-18. All we need
				 // to change is hence to perform the
				 // cell-related operations only on
				 // the process with the correct
				 // ID. The second difference is the
				 // way we calculate the maximum
				 // value. Before, we could simply
				 // have a <code>double</code>
				 // variable that we checked against
				 // on each quadrature point for each
				 // cell. Now, we have to be a bit
				 // more careful since each processor
				 // only operates on a subset of
				 // cells. What we do is to first let
				 // each processor calculate the
				 // maximum among its cells, and then
				 // do a global communication
				 // operation called
				 // <code>MaxAll</code> that searches
				 // for the maximum value among all
				 // the maximum values of the
				 // individual processors. MPI
				 // provides such a call, but it's
				 // even simpler to use the respective
				 // function of the MPI
				 // communicator object since that
				 // will do the right thing even if we
				 // work without MPI and on a single
				 // machine only. The call to
				 // <code>MaxAll</code> needs three
				 // arguments, namely the local
				 // maximum (input), a field for the
				 // global maximum (output), and an
				 // integer value one that says that
				 // we only work on one double.
template <int dim>
double BoussinesqFlowProblem<dim>::get_maximal_velocity () const
{
  const QIterated<dim> quadrature_formula (QTrapez<1>(),
					   stokes_degree+1);
  const unsigned int n_q_points = quadrature_formula.size();

  FEValues<dim> fe_values (stokes_fe, quadrature_formula, update_values);
  std::vector<Tensor<1,dim> > velocity_values(n_q_points);

  const FEValuesExtractors::Vector velocities (0);

  double max_local_velocity = 0;

  typename DoFHandler<dim>::active_cell_iterator
    cell = stokes_dof_handler.begin_active(),
    endc = stokes_dof_handler.end();
  for (; cell!=endc; ++cell)
    if (cell->subdomain_id() ==
	Utilities::System::get_this_mpi_process(MPI_COMM_WORLD))
      {
	fe_values.reinit (cell);
	fe_values[velocities].get_function_values (stokes_solution,
						   velocity_values);

	for (unsigned int q=0; q<n_q_points; ++q)
	  max_local_velocity = std::max (max_local_velocity,
					 velocity_values[q].norm());
      }

  double max_velocity = 0.;
#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
  MPI_Allreduce (&max_local_velocity, &max_velocity, 1, MPI_DOUBLE,
		 MPI_MAX, MPI_COMM_WORLD);
#else
  max_velocity = max_local_velocity;
#endif

  return max_velocity;
}



				 // Again, this is only a slightly
				 // modified version of the respective
				 // function in step-31. What is new is
				 // that each processor works on its
				 // partition of cells, and gets a minimum
				 // and maximum temperature on that
				 // partition. Two global communication
				 // steps synchronize the data among the
				 // processors.
template <int dim>
std::pair<double,double>
BoussinesqFlowProblem<dim>::get_extrapolated_temperature_range () const
{
  const QIterated<dim> quadrature_formula (QTrapez<1>(),
					   temperature_degree);
  const unsigned int n_q_points = quadrature_formula.size();

  FEValues<dim> fe_values (temperature_fe, quadrature_formula,
                           update_values);
  std::vector<double> old_temperature_values(n_q_points);
  std::vector<double> old_old_temperature_values(n_q_points);

				   // This presets the minimum with a bigger
				   // and the maximum with a smaller number
				   // than one that is going to appear. Will
				   // be overwritten in the cell loop or in
				   // the communication step at the
				   // latest.
  double min_local_temperature = std::numeric_limits<double>::max(),
	 max_local_temperature = -std::numeric_limits<double>::max();

  if (timestep_number != 0)
    {
      typename DoFHandler<dim>::active_cell_iterator
	cell = temperature_dof_handler.begin_active(),
	endc = temperature_dof_handler.end();
      for (; cell!=endc; ++cell)
	if (cell->subdomain_id() ==
	    Utilities::System::get_this_mpi_process(MPI_COMM_WORLD))
	  {
	    fe_values.reinit (cell);
	    fe_values.get_function_values (old_temperature_solution,
					   old_temperature_values);
	    fe_values.get_function_values (old_old_temperature_solution,
					   old_old_temperature_values);

	    for (unsigned int q=0; q<n_q_points; ++q)
	      {
		const double temperature =
		  (1. + time_step/old_time_step) * old_temperature_values[q]-
		  time_step/old_time_step * old_old_temperature_values[q];

		min_local_temperature = std::min (min_local_temperature,
						  temperature);
		max_local_temperature = std::max (max_local_temperature,
						  temperature);
	      }
	  }
    }
  else
    {
      typename DoFHandler<dim>::active_cell_iterator
	cell = temperature_dof_handler.begin_active(),
	endc = temperature_dof_handler.end();
      for (; cell!=endc; ++cell)
	if (cell->subdomain_id() ==
	    Utilities::System::get_this_mpi_process(MPI_COMM_WORLD))
	  {
	    fe_values.reinit (cell);
	    fe_values.get_function_values (old_temperature_solution,
					   old_temperature_values);

	    for (unsigned int q=0; q<n_q_points; ++q)
	      {
		const double temperature = old_temperature_values[q];

		min_local_temperature = std::min (min_local_temperature,
						  temperature);
		max_local_temperature = std::max (max_local_temperature,
						  temperature);
	      }
	  }
    }
  
  double min_temperature, max_temperature;
#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
  MPI_Allreduce (&max_local_temperature, &max_temperature, 1, MPI_DOUBLE,
		 MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce (&min_local_temperature, &min_temperature, 1, MPI_DOUBLE,
		 MPI_MIN, MPI_COMM_WORLD);
#else
  min_temperature = min_local_temperature;
  max_temperature = max_local_temperature;
#endif

  return std::make_pair(min_temperature, max_temperature);
}



				 // The function that calculates the
				 // viscosity is purely local, so this is
				 // the same code as in step-31.
template <int dim>
double
BoussinesqFlowProblem<dim>::
compute_viscosity (const std::vector<double>          &old_temperature,
		   const std::vector<double>          &old_old_temperature,
		   const std::vector<Tensor<1,dim> >  &old_temperature_grads,
		   const std::vector<Tensor<1,dim> >  &old_old_temperature_grads,
		   const std::vector<double>          &old_temperature_laplacians,
		   const std::vector<double>          &old_old_temperature_laplacians,
		   const std::vector<Tensor<1,dim> >  &old_velocity_values,
		   const std::vector<Tensor<1,dim> >  &old_old_velocity_values,
		   const std::vector<SymmetricTensor<2,dim> >  &old_strain_rates,
		   const std::vector<SymmetricTensor<2,dim> >  &old_old_strain_rates,
		   const double                        global_u_infty,
		   const double                        global_T_variation,
		   const double                        cell_diameter) const
{
  const double beta = 0.026 * dim;
  const double alpha = 1;

  if (global_u_infty == 0)
    return 5e-3 * cell_diameter;

  const unsigned int n_q_points = old_temperature.size();

  double max_residual = 0;
  double max_velocity = 0;

  for (unsigned int q=0; q < n_q_points; ++q)
    {
      const Tensor<1,dim> u = (old_velocity_values[q] +
			       old_old_velocity_values[q]) / 2;

      const SymmetricTensor<2,dim> strain_rate = (old_strain_rates[q] +
						  old_old_strain_rates[q]) / 2;
      
      const double T = (old_temperature[q] + old_old_temperature[q]) / 2;
      const double dT_dt = (old_temperature[q] - old_old_temperature[q])
			   / old_time_step;
      const double u_grad_T = u * (old_temperature_grads[q] +
				   old_old_temperature_grads[q]) / 2;

      const double kappa_Delta_T = EquationData::kappa
				   * (old_temperature_laplacians[q] +
				      old_old_temperature_laplacians[q]) / 2;
      const double gamma
	= ((EquationData::radiogenic_heating * EquationData::density(T)
	    +
	    2 * EquationData::eta * strain_rate * strain_rate) /
	   (EquationData::density(T) * EquationData::specific_heat));      

      const double residual
	= std::abs((dT_dt + u_grad_T - kappa_Delta_T - gamma) *
		   std::pow((old_temperature[q]+old_old_temperature[q]) / 2,
			    alpha-1.));

      max_residual = std::max (residual,        max_residual);
      max_velocity = std::max (std::sqrt (u*u), max_velocity);
    }

  if (timestep_number == 0)
    return beta * max_velocity * cell_diameter;
  else
    {
      Assert (old_time_step > 0, ExcInternalError());

      const double c_R = 0.11;
      const double global_scaling = c_R * global_u_infty * global_T_variation *
				    std::pow(global_Omega_diameter, alpha - 2.);

      return (beta *
	      max_velocity *
	      std::min (cell_diameter,
			std::pow(cell_diameter,alpha) * max_residual /
			global_scaling));
    }
}



				 // This function is new compared to
				 // step-31. What is does is to re-implement
				 // the library function
				 // <code>VectorTools::project()</code> for
				 // an MPI-based parallelization, a function
				 // we used for generating an initial vector
				 // for temperature based on some initial
				 // function. The library function only
				 // works with shared memory but doesn't
				 // know how to utilize multiple machines
				 // coupled through MPI to compute the
				 // projected solution. If run with
				 // more than one MPI process, this would
				 // mean that each processor projects the
				 // whole field, which is clearly not very
				 // efficient. The details of a
				 // <code>project()</code> function are not
				 // very difficult. All we do is to use a
				 // mass matrix and put the evaluation of
				 // the initial value function on the right
				 // hand side. The mass matrix for
				 // temperature we can simply generate using
				 // the respective assembly function, so all
				 // we need to do here is to create the
				 // right hand side and do a CG solve. The
				 // assembly function does a loop over all
				 // cells and evaluates the function in the
				 // <code>EquationData</code> namespace, and
				 // does this only on cells pertaining to
				 // the respective processor. The
				 // implementation of this assembly differs
				 // from the assembly we do for the
				 // principal assembly functions further
				 // down (which include thread-based
				 // parallelization with the WorkStream
				 // concept). Here we chose to keep things
				 // simple (keeping in mind that this function
				 // is also only called at the beginning of
				 // the program, not every time step), and
				 // generating that right hand
				 // side is cheap anyway so we won't even
				 // notice that this part is not parallized
				 // by threads.
				 //
				 // Regarding the implementation of
				 // inhomogeneous Dirichlet boundary
				 // conditions: Since we use the temperature
				 // ConstraintMatrix, we can apply the
				 // boundary conditions directly when
				 // building the respective matrix and right
				 // hand side. In this case, the boundary
				 // conditions are inhomogeneous, which
				 // makes this procedure somewhat
				 // tricky. Remember that we get the matrix
				 // from some other function. However, the
				 // correct imposition of boundary
				 // conditions needs the matrix data we work
				 // on plus the right hand side
				 // simultaneously, since the right hand
				 // side is created by Gaussian elimination
				 // on the matrix rows. In order to not
				 // introduce the matrix assembly at this
				 // place, but still having the matrix data
				 // available, we choose to create a dummy
				 // matrix <code>matrix_for_bc</code> that
				 // we only fill with data when we need it
				 // for imposing boundary conditions. These
				 // positions are exactly those where we
				 // have an inhomogeneous entry in the
				 // ConstraintMatrix. There are only a few
				 // such positions (on the boundary dofs),
				 // so it is still much cheaper to use this
				 // function than to create the full matrix
				 // here. To implement this, we ask the
				 // constraint matrix whether the dof under
				 // consideration is inhomogeneously
				 // constraint. In that case, we generate
				 // the respective matrix column that we
				 // need for creating the correct right hand
				 // side. Note that this (manually
				 // generated) matrix entry needs to be
				 // exactly the entry that we would fill the
				 // matrix with &mdash; otherwise, this will
				 // not work.
template <int dim>
void BoussinesqFlowProblem<dim>::project_temperature_field ()
{
  assemble_temperature_matrix ();

  QGauss<dim> quadrature(temperature_degree+2);
  UpdateFlags update_flags = UpdateFlags(update_values   |
					 update_quadrature_points |
					 update_JxW_values);
  FEValues<dim> fe_values (temperature_fe, quadrature, update_flags);

  const unsigned int dofs_per_cell = fe_values.dofs_per_cell,
		     n_q_points    = fe_values.n_quadrature_points;

  std::vector<unsigned int> local_dof_indices (dofs_per_cell);
  Vector<double> cell_vector (dofs_per_cell);
  FullMatrix<double> matrix_for_bc (dofs_per_cell, dofs_per_cell);

  typename DoFHandler<dim>::active_cell_iterator
    cell = temperature_dof_handler.begin_active(),
    endc = temperature_dof_handler.end();

  std::vector<double> rhs_values(n_q_points);

  TrilinosWrappers::MPI::Vector
    rhs (temperature_mass_matrix.row_partitioner()),
    solution (temperature_mass_matrix.row_partitioner());

  for (; cell!=endc; ++cell)
    if (cell->subdomain_id() ==
	Utilities::System::get_this_mpi_process(MPI_COMM_WORLD))
      {
	cell->get_dof_indices (local_dof_indices);
	fe_values.reinit(cell);

	EquationData::TemperatureInitialValues<dim>().value_list
	  (fe_values.get_quadrature_points(), rhs_values);

	cell_vector = 0;
	matrix_for_bc = 0;
	for (unsigned int point=0; point<n_q_points; ++point)
	  for (unsigned int i=0; i<dofs_per_cell; ++i)
	    {
	      cell_vector(i) += rhs_values[point] *
		                fe_values.shape_value(i,point) *
		                fe_values.JxW(point);
	      if (temperature_constraints.is_inhomogeneously_constrained(local_dof_indices[i]))
		{
		  for (unsigned int j=0; j<dofs_per_cell; ++j)
		    matrix_for_bc(j,i) += fe_values.shape_value(i,point) *
					  fe_values.shape_value(j,point) *
		                          fe_values.JxW(point);
		}
	    }

	temperature_constraints.distribute_local_to_global (cell_vector,
							    local_dof_indices,
							    rhs,
							    matrix_for_bc);
      }

  rhs.compress (Add);

  SolverControl solver_control(5*rhs.size(), 1e-12*rhs.l2_norm());
  SolverCG<TrilinosWrappers::MPI::Vector> cg(solver_control);

  TrilinosWrappers::PreconditionJacobi preconditioner_mass;
  preconditioner_mass.initialize(temperature_mass_matrix, 1.3);

  cg.solve (temperature_mass_matrix, solution, rhs, preconditioner_mass);

  temperature_constraints.distribute (solution);

//  old_temperature_solution = solution;
  old_temperature_solution.reinit(solution, false, true);
}




				 // @sect4{The BoussinesqFlowProblem setup functions}

				 // The following three functions set
				 // up the Stokes matrix, the matrix
				 // used for the Stokes
				 // preconditioner, and the
				 // temperature matrix. The code is
				 // mostly the same as in step-31, but
				 // it has been broken out into three
				 // functions of their own for
				 // simplicity, but also so that they
				 // can easily be run in %parallel on
				 // multiple threads (unless we are
				 // running with MPI, in which case
				 // this is not possible, as explained
				 // in the introduction).
				 //
				 // The main functional difference
				 // between the code here and that in
				 // step-31 is that the matrices we
				 // want to set up are distributed
				 // across multiple processors. Since
				 // we still want to build up the
				 // sparsity pattern first for
				 // efficiency reasons, we could
				 // continue to build the
				 // <i>entire</i> sparsity pattern as
				 // a
				 // BlockCompressedSimpleSparsityPattern,
				 // as we did in step-31. However,
				 // that would be inefficient: every
				 // processor would build the same
				 // sparsity pattern, but only
				 // initialize a small part of the
				 // matrix using it.
				 //
				 // Rather, we use an object of type
				 // TrilinosWrappers::BlockSparsityPattern,
				 // which is (obviously) a wrapper
				 // around a sparsity pattern object
				 // provided by Trilinos. The
				 // advantage is that the Trilinos
				 // sparsity pattern class can
				 // communicate across multiple
				 // processors: if this processor
				 // fills in all the nonzero entries
				 // that result from the cells it
				 // owns, and every other processor
				 // does so as well, then at the end
				 // after some MPI communication
				 // initiated by the
				 // <code>compress()</code> call, we
				 // will have the globally assembled
				 // sparsity pattern available with
				 // which the global matrix can be
				 // initialized.
				 //
				 // The only other change we need to
				 // make is to tell the
				 // DoFTools::make_sparsity_pattern
				 // function that it is only supposed
				 // to work on a subset of cells,
				 // namely the ones whose
				 // <code>subdomain_id</code> equals
				 // the number of the current
				 // processor, and to ignore all other
				 // cells.
				 //
				 // This strategy is replicated across
				 // all three of the following
				 // functions.
				 //
				 // Note that Trilinos matrices store the
 				 // information contained in the sparsity
 				 // patterns, so we can safely release the
 				 // <code>sp</code> variable once the matrix
 				 // has been given the sparsity structure.
template <int dim>
void BoussinesqFlowProblem<dim>::
  setup_stokes_matrix (const std::vector<IndexSet> &stokes_partitioning)
{
  stokes_matrix.clear ();

  TrilinosWrappers::BlockSparsityPattern sp (stokes_partitioning,
					     MPI_COMM_WORLD);

  Table<2,DoFTools::Coupling> coupling (dim+1, dim+1);

  for (unsigned int c=0; c<dim+1; ++c)
    for (unsigned int d=0; d<dim+1; ++d)
      if (! ((c==dim) && (d==dim)))
	coupling[c][d] = DoFTools::always;
      else
	coupling[c][d] = DoFTools::none;

  DoFTools::make_sparsity_pattern (stokes_dof_handler,
				   coupling, sp,
				   stokes_constraints, false,
				   Utilities::System::
				   get_this_mpi_process(MPI_COMM_WORLD));
  sp.compress();

  stokes_matrix.reinit (sp);
}



template <int dim>
void BoussinesqFlowProblem<dim>::
  setup_stokes_preconditioner (const std::vector<IndexSet> &stokes_partitioning)
{
  Amg_preconditioner.reset ();
  Mp_preconditioner.reset ();

  stokes_preconditioner_matrix.clear ();

  TrilinosWrappers::BlockSparsityPattern sp (stokes_partitioning,
					     MPI_COMM_WORLD);

  Table<2,DoFTools::Coupling> coupling (dim+1, dim+1);
  for (unsigned int c=0; c<dim+1; ++c)
    for (unsigned int d=0; d<dim+1; ++d)
      if (c == d)
	coupling[c][d] = DoFTools::always;
      else
	coupling[c][d] = DoFTools::none;

  DoFTools::make_sparsity_pattern (stokes_dof_handler,
				   coupling, sp,
				   stokes_constraints, false,
				   Utilities::System::
				   get_this_mpi_process(MPI_COMM_WORLD));
  sp.compress();

  stokes_preconditioner_matrix.reinit (sp);
}


template <int dim>
void BoussinesqFlowProblem<dim>::
  setup_temperature_matrices (const IndexSet &temperature_partitioner)
{
  T_preconditioner.reset ();
  temperature_mass_matrix.clear ();
  temperature_stiffness_matrix.clear ();
  temperature_matrix.clear ();

  TrilinosWrappers::SparsityPattern sp (temperature_partitioner,
					MPI_COMM_WORLD);
  DoFTools::make_sparsity_pattern (temperature_dof_handler, sp,
				   temperature_constraints, false,
				   Utilities::System::
				   get_this_mpi_process(MPI_COMM_WORLD));
  sp.compress();

  temperature_matrix.reinit (sp);
  temperature_mass_matrix.reinit (sp);
  temperature_stiffness_matrix.reinit (sp);
}



				 // The remainder of the setup function
				 // (after splitting out the three functions
				 // above) mostly has to deal with the
				 // things we need to do for parallelization
				 // across processors. In particular, at the
				 // top it calls
				 // GridTools::partition_triangulation to
				 // subdivide all cells into subdomains of
				 // roughly equal size and roughly minimal
				 // interface length (using METIS). We then
				 // distribute degrees of freedom for Stokes
				 // and temperature DoFHandler objects, and
				 // re-sort them in such a way that all
				 // degrees of freedom associated with
				 // subdomain zero come before all those
				 // associated with subdomain one, etc. For
				 // the Stokes part, this entails, however,
				 // that velocities and pressures become
				 // intermixed, but this is trivially solved
				 // by sorting again by blocks; it is worth
				 // noting that this latter operation leaves
				 // the relative ordering of all velocities
				 // and pressures alone, i.e. within the
				 // velocity block we will still have all
				 // those associated with subdomain zero
				 // before all velocities associated with
				 // subdomain one, etc. This is important
				 // since we store each of the blocks of
				 // this matrix distributed across all
				 // processors and want this to be done in
				 // such a way that each processor stores
				 // that part of the matrix that is roughly
				 // equal to the degrees of freedom located
				 // on those cells that it will actually
				 // work on. Note how we set boundary
				 // conditions on the temperature by using
				 // the ConstraintMatrix object.
				 //
				 // After this, we have to set up the
				 // various partitioners (of type
				 // <code>IndexSet</code>, see the
				 // introduction) that describe which parts
				 // of each matrix or vector will be stored
				 // where, then call the functions that
				 // actually set up the matrices
				 // (concurrently if not using MPI
				 // but sequentially otherwise, as explained
				 // in the introduction), and at the end also
				 // resize the various vectors we keep
				 // around in this program. We given those
				 // vectors the correct size using the
				 // aforementioned Epetra_Map. Most of the
				 // vectors are actually localized, i.e.,
				 // they store all dofs in the problem on
				 // each processor. In that case, the only
				 // information that is used is the global
				 // size. This is different for the two
				 // right hand side vectors, which are
				 // distributed ones, see also the class
				 // declaration.
				 //
				 // Note how this function enters and leaves
				 // a timed section so that we can get a
				 // time report at the end of the
				 // program. Note also the use of the
				 // <code>pcout</code> variable: to every
				 // process it looks like we can write to
				 // screen, but only the output of the first
				 // processor actually ends up somewhere. We
				 // could of course have achieved the same
				 // effect by writing to
				 // <code>std::cout</code> but would then
				 // have had to guard every access to that
				 // stream by something like <code>if
				 // (Utilities:: System::
				 // get_this_mpi_process
				 // (MPI_COMM_WORLD) == 0)</code>,
				 // hardly a pretty solution.
template <int dim>
void BoussinesqFlowProblem<dim>::setup_dofs ()
{
  computing_timer.enter_section("Setup dof systems");

  std::vector<unsigned int> stokes_sub_blocks (dim+1,0);
  stokes_sub_blocks[dim] = 1;
  stokes_dof_handler.distribute_dofs (stokes_fe);
  DoFRenumbering::component_wise (stokes_dof_handler, stokes_sub_blocks);

  temperature_dof_handler.distribute_dofs (temperature_fe);

  std::vector<unsigned int> stokes_dofs_per_block (2);
  DoFTools::count_dofs_per_block (stokes_dof_handler, stokes_dofs_per_block,
				  stokes_sub_blocks);

  const unsigned int n_u = stokes_dofs_per_block[0],
                     n_p = stokes_dofs_per_block[1],
		     n_T = temperature_dof_handler.n_dofs();

  pcout << "Number of active cells: "
	<< triangulation.n_global_active_cells()
	<< " (on "
	<< triangulation.n_levels()
	<< " levels)"
	<< std::endl
	<< "Number of degrees of freedom: "
	<< n_u + n_p + n_T
	<< " (" << n_u << '+' << n_p << '+'<< n_T <<')'
	<< std::endl
	<< std::endl;



  std::vector<IndexSet> stokes_partitioning, stokes_relevant_partitioning;
  IndexSet temperature_partitioning (n_T), temperature_relevant_partitioning (n_T);
  IndexSet stokes_relevant_set;
  {
    const unsigned int my_id =
      Utilities::System::get_this_mpi_process(MPI_COMM_WORLD);
    IndexSet stokes_index_set = stokes_dof_handler.locally_owned_dofs();
    stokes_partitioning.push_back(stokes_index_set.get_view(0,n_u));
    stokes_partitioning.push_back(stokes_index_set.get_view(n_u,n_u+n_p));

    DoFTools::extract_locally_relevant_dofs (stokes_dof_handler,
					     stokes_relevant_set);
    stokes_relevant_partitioning.push_back(stokes_relevant_set.get_view(0,n_u));
    stokes_relevant_partitioning.push_back(stokes_relevant_set.get_view(n_u,n_u+n_p));

    temperature_partitioning = temperature_dof_handler.locally_owned_dofs();
    DoFTools::extract_locally_relevant_dofs (temperature_dof_handler,
					     temperature_relevant_partitioning);
  }

  {

    stokes_constraints.clear ();
//    IndexSet stokes_la;
//    DoFTools::extract_locally_active_dofs (stokes_dof_handler,
//					   stokes_la);
    stokes_constraints.reinit(stokes_relevant_set);

    DoFTools::make_hanging_node_constraints (stokes_dof_handler,
					     stokes_constraints);

    std::vector<bool> velocity_mask (dim+1, true);
    velocity_mask[dim] = false;
    VectorTools::interpolate_boundary_values (stokes_dof_handler,
					      0,
					      ZeroFunction<dim>(dim+1),
					      stokes_constraints,
					      velocity_mask);

    std::set<unsigned char> no_normal_flux_boundaries;
    no_normal_flux_boundaries.insert (1);
    VectorTools::compute_no_normal_flux_constraints (stokes_dof_handler, 0,
						     no_normal_flux_boundaries,
						     stokes_constraints);
    stokes_constraints.close ();
  }
  {
    temperature_constraints.clear ();
    temperature_constraints.reinit(temperature_relevant_partitioning);//temp_locally_active);

    DoFTools::make_hanging_node_constraints (temperature_dof_handler,
					     temperature_constraints);
    VectorTools::interpolate_boundary_values (temperature_dof_handler,
					      0,
					      EquationData::TemperatureInitialValues<dim>(),
					      temperature_constraints);
    VectorTools::interpolate_boundary_values (temperature_dof_handler,
					      1,
					      EquationData::TemperatureInitialValues<dim>(),
					      temperature_constraints);
    temperature_constraints.close ();
  }

  if (Utilities::System::job_supports_mpi() == false)
    {
      Threads::TaskGroup<> tasks;
      tasks += Threads::new_task (&BoussinesqFlowProblem<dim>::setup_stokes_matrix,
				  *this,
				  stokes_partitioning);
      tasks += Threads::new_task (&BoussinesqFlowProblem<dim>::setup_stokes_preconditioner,
				  *this,
				  stokes_partitioning);
      tasks += Threads::new_task (&BoussinesqFlowProblem<dim>::setup_temperature_matrices,
				  *this,
				  temperature_partitioning);
      tasks.join_all ();
    }
  else
    {
      setup_stokes_matrix (stokes_partitioning);
      setup_stokes_preconditioner (stokes_partitioning);
      setup_temperature_matrices (temperature_partitioning);
    }

  stokes_rhs.reinit (stokes_partitioning, MPI_COMM_WORLD);
  stokes_solution.reinit (stokes_relevant_partitioning, MPI_COMM_WORLD);
  old_stokes_solution.reinit (stokes_solution);

  temperature_rhs.reinit (temperature_partitioning, MPI_COMM_WORLD);
  temperature_solution.reinit (temperature_relevant_partitioning, MPI_COMM_WORLD);
  old_temperature_solution.reinit (temperature_solution);
  old_old_temperature_solution.reinit (temperature_solution);

  rebuild_stokes_matrix              = true;
  rebuild_stokes_preconditioner      = true;
  rebuild_temperature_matrices       = true;
  rebuild_temperature_preconditioner = true;

  computing_timer.exit_section();
}



				 // @sect4{The BoussinesqFlowProblem assembly functions}
				 //
				 // Following the discussion in the
				 // introduction and in the @ref threads
				 // module, we split the assembly functions
				 // into different parts:
				 //
				 // <ul>
				 // <li> The local calculations of matrices
				 // and right hand sides, given a certain cell
				 // as input (these functions are named
				 // <code>local_assemble_*</code> below). The
				 // resulting function is, in other words,
				 // essentially the body of the loop over all
				 // cells in step-31. Note, however, that
				 // these functions store the result from the
				 // local calculations in variables of classes
				 // from the CopyData namespace.
				 //
				 // <li>These objects are then given to the
				 // second step which writes the local data
				 // into the global data structures (these
				 // functions are named
				 // <code>copy_local_to_global_*</code>
				 // below). These functions are pretty
				 // trivial.
				 //
				 // <li>These two subfunctions are then used
				 // in the respective assembly routine (called
				 // <code>assemble_*</code> below), where a
				 // WorkStream object is set up and runs over
				 // all the cells that belong to the
				 // processor's subdomain.
				 // </ul>

				 // @sect5{Stokes preconditioner assembly}
				 //
				 // Let us start with the functions that
				 // builds the Stokes preconditioner. The
				 // first two of these are pretty trivial,
				 // given the discussion above. Note in
				 // particular that the main point in using
				 // the scratch data object is that we want to
				 // avoid allocating any objects on the free
				 // space each time we visit a new cell. As a
				 // consequence, the assembly function below
				 // only has automatic local variables, and
				 // everything else is accessed through the
				 // scratch data object, which is allocated
				 // only once before we start the loop over
				 // all cells:
template <int dim>
void
BoussinesqFlowProblem<dim>::
local_assemble_stokes_preconditioner (const typename DoFHandler<dim>::active_cell_iterator &cell,
				      Assembly::Scratch::StokesPreconditioner<dim> &scratch,
				      Assembly::CopyData::StokesPreconditioner<dim> &data)
{
  const unsigned int   dofs_per_cell   = stokes_fe.dofs_per_cell;
  const unsigned int   n_q_points      = scratch.stokes_fe_values.n_quadrature_points;

  const FEValuesExtractors::Vector velocities (0);
  const FEValuesExtractors::Scalar pressure (dim);

  scratch.stokes_fe_values.reinit (cell);
  cell->get_dof_indices (data.local_dof_indices);

  data.local_matrix = 0;

  for (unsigned int q=0; q<n_q_points; ++q)
    {
      for (unsigned int k=0; k<dofs_per_cell; ++k)
	{
	  scratch.grad_phi_u[k] = scratch.stokes_fe_values[velocities].gradient(k,q);
	  scratch.phi_p[k]      = scratch.stokes_fe_values[pressure].value (k, q);
	}

      for (unsigned int i=0; i<dofs_per_cell; ++i)
	for (unsigned int j=0; j<dofs_per_cell; ++j)
	  data.local_matrix(i,j) += (EquationData::eta *
				     scalar_product (scratch.grad_phi_u[i],
						     scratch.grad_phi_u[j])
				     +
				     (1./EquationData::eta) *
				     EquationData::pressure_scaling *
				     EquationData::pressure_scaling *
				     (scratch.phi_p[i] * scratch.phi_p[j]))
				    * scratch.stokes_fe_values.JxW(q);
    }
}



template <int dim>
void
BoussinesqFlowProblem<dim>::
copy_local_to_global_stokes_preconditioner (const Assembly::CopyData::StokesPreconditioner<dim> &data)
{
  stokes_constraints.distribute_local_to_global (data.local_matrix,
						 data.local_dof_indices,
						 stokes_preconditioner_matrix);
}



				 // When we create the WorkStream, we modify
				 // the start and end iterator into a
				 // so-called <code>SubdomainFilter</code>
				 // that tells the individual processes which
				 // cells to work on. This is exactly the case
				 // discussed in the introduction. Note how we
				 // use the construct
				 // <code>std_cxx1x::bind</code> to create a
				 // function object that is compatible with
				 // the WorkStream class. It uses placeholders
				 // <code>_1, _2, _3</code> for the local
				 // assembly function that specify cell,
				 // scratch data, and copy data, as well as
				 // the placeholder <code>_1</code> for the
				 // copy function that expects the data to be
				 // written into the global matrix. On the
				 // other hand, the implicit zeroth argument
				 // of member functions (namely the
				 // <code>this</code> pointer of the object on
				 // which that member function is to operate
				 // on) is <i>bound</i> to the
				 // <code>this</code> pointer of the current
				 // function. The WorkStream class, as a
				 // consequence, does not need to know
				 // anything about the object these functions
				 // work on.
				 //
				 // When the
				 // WorkStream is executed, it will create
				 // several local assembly routines of the
				 // first kind for several cells and let
				 // some available processors work on
				 // them. The function that needs to be
				 // synchronized, i.e., the write operation
				 // into the global matrix, however, is
				 // executed by only one thread at a time in
				 // the prescribed order. Of course, this
				 // only holds for the parallelization on a
				 // single MPI process. Different MPI
				 // processes will have their own WorkStream
				 // objects and do that work completely
				 // independently. In a distributed
				 // calculation, some data will accumulate
				 // at degrees of freedom that are not owned
				 // by the respective processor. It would be
				 // inefficient to send data around every
				 // time we encounter such a dof. What
				 // happens instead is that the Trilinos
				 // sparse matrix will keep that data and
				 // send it to the owner at the end of
				 // assembly, by calling the
				 // <code>compress()</code> command.
template <int dim>
void
BoussinesqFlowProblem<dim>::assemble_stokes_preconditioner ()
{
  stokes_preconditioner_matrix = 0;

  const QGauss<dim> quadrature_formula(stokes_degree+2);

  typedef
    FilteredIterator<typename DoFHandler<dim>::active_cell_iterator>
    SubdomainFilter;

  WorkStream::
    run (SubdomainFilter (IteratorFilters::SubdomainEqualTo
			  (Utilities::System::get_this_mpi_process(MPI_COMM_WORLD)),
			  stokes_dof_handler.begin_active()),
	 SubdomainFilter (IteratorFilters::SubdomainEqualTo
			  (Utilities::System::get_this_mpi_process(MPI_COMM_WORLD)),
			  stokes_dof_handler.end()),
	 std_cxx1x::bind (&BoussinesqFlowProblem<dim>::
			  local_assemble_stokes_preconditioner,
			  this,
			  _1,
			  _2,
			  _3),
	 std_cxx1x::bind (&BoussinesqFlowProblem<dim>::
			  copy_local_to_global_stokes_preconditioner,
			  this,
			  _1),
	 Assembly::Scratch::
	 StokesPreconditioner<dim> (stokes_fe, quadrature_formula,
				    update_JxW_values |
				    update_values |
				    update_gradients),
	 Assembly::CopyData::
	 StokesPreconditioner<dim> (stokes_fe));

  stokes_preconditioner_matrix.compress();
}



				 // The final function in this block initiates
				 // assemble of the Stokes preconditioner
				 // matrix and then builds the Stokes
				 // preconditioner. It is mostly the same as
				 // in the serial case. The only difference to
				 // step-31 is that we use an ILU
				 // preconditioner for the pressure mass
				 // matrix instead of IC, as discussed in the
				 // introduction.
template <int dim>
void
BoussinesqFlowProblem<dim>::build_stokes_preconditioner ()
{
  if (rebuild_stokes_preconditioner == false)
    return;

  computing_timer.enter_section ("   Build Stokes preconditioner");
  pcout << "   Rebuilding Stokes preconditioner..." << std::flush;

  assemble_stokes_preconditioner ();

  std::vector<std::vector<bool> > constant_modes;
  std::vector<bool>  velocity_components (dim+1,true);
  velocity_components[dim] = false;
  DoFTools::extract_constant_modes (stokes_dof_handler, velocity_components,
				    constant_modes);

  Mp_preconditioner  = std_cxx1x::shared_ptr<TrilinosWrappers::PreconditionILU>
		       (new TrilinosWrappers::PreconditionILU());
  Amg_preconditioner = std_cxx1x::shared_ptr<TrilinosWrappers::PreconditionAMG>
		       (new TrilinosWrappers::PreconditionAMG());

  TrilinosWrappers::PreconditionAMG::AdditionalData Amg_data;
  Amg_data.constant_modes = constant_modes;
  Amg_data.elliptic = true;
  Amg_data.higher_order_elements = true;
  Amg_data.smoother_sweeps = 2;
//  Amg_data.aggregation_threshold = 0.02;

  Mp_preconditioner->initialize (stokes_preconditioner_matrix.block(1,1));
  Amg_preconditioner->initialize (stokes_preconditioner_matrix.block(0,0),
				  Amg_data);

  rebuild_stokes_preconditioner = false;

 pcout << std::endl;
 computing_timer.exit_section();
}

				 // @sect5{Stokes system assembly}

				 // The next three functions implement the
				 // assembly of the Stokes system, again
				 // split up into a part performing local
				 // calculations, one for writing the local
				 // data into the global matrix and vector,
				 // and one for actually running the loop
				 // over all cells with the help of the
				 // WorkStream class. Note that the assembly
				 // of the Stokes matrix needs only to be
				 // done in case we have changed the
				 // mesh. Otherwise, just the
				 // (temperature-dependent) right hand side
				 // needs to be calculated here. Since we
				 // are working with distributed matrices
				 // and vectors, we have to call the
				 // respective <code>compress()</code>
				 // functions in the end of the assembly in
				 // order to send non-local data to the
				 // owner process.
template <int dim>
void
BoussinesqFlowProblem<dim>::
local_assemble_stokes_system (const typename DoFHandler<dim>::active_cell_iterator &cell,
			      Assembly::Scratch::StokesSystem<dim> &scratch,
			      Assembly::CopyData::StokesSystem<dim> &data)
{
  const unsigned int dofs_per_cell = scratch.stokes_fe_values.get_fe().dofs_per_cell;
  const unsigned int n_q_points    = scratch.stokes_fe_values.n_quadrature_points;

  const FEValuesExtractors::Vector velocities (0);
  const FEValuesExtractors::Scalar pressure (dim);

  scratch.stokes_fe_values.reinit (cell);

  typename DoFHandler<dim>::active_cell_iterator
    temperature_cell (&triangulation,
		      cell->level(),
		      cell->index(),
		      &temperature_dof_handler);
  scratch.temperature_fe_values.reinit (temperature_cell);

  if (rebuild_stokes_matrix)
    data.local_matrix = 0;
  data.local_rhs = 0;

  scratch.temperature_fe_values.get_function_values (old_temperature_solution,
						     scratch.old_temperature_values);

  for (unsigned int q=0; q<n_q_points; ++q)
    {
      const double old_temperature = scratch.old_temperature_values[q];

      for (unsigned int k=0; k<dofs_per_cell; ++k)
	{
	  scratch.phi_u[k] = scratch.stokes_fe_values[velocities].value (k,q);
	  if (rebuild_stokes_matrix)
	    {
	      scratch.grads_phi_u[k] = scratch.stokes_fe_values[velocities].symmetric_gradient(k,q);
	      scratch.div_phi_u[k]   = scratch.stokes_fe_values[velocities].divergence (k, q);
	      scratch.phi_p[k]       = scratch.stokes_fe_values[pressure].value (k, q);
	    }
	}

      if (rebuild_stokes_matrix)
	for (unsigned int i=0; i<dofs_per_cell; ++i)
	  for (unsigned int j=0; j<dofs_per_cell; ++j)
	    data.local_matrix(i,j) += (EquationData::eta * 2 *
				       (scratch.grads_phi_u[i] * scratch.grads_phi_u[j])
				       - (EquationData::pressure_scaling *
					  scratch.div_phi_u[i] * scratch.phi_p[j])
				       - (EquationData::pressure_scaling *
					  scratch.phi_p[i] * scratch.div_phi_u[j]))
				      * scratch.stokes_fe_values.JxW(q);

      const Tensor<1,dim>
	gravity = EquationData::gravity_vector (scratch.stokes_fe_values
						.quadrature_point(q));

      for (unsigned int i=0; i<dofs_per_cell; ++i)
	data.local_rhs(i) += (EquationData::density(old_temperature) *
			      gravity  *
			      scratch.phi_u[i]) *
			     scratch.stokes_fe_values.JxW(q);
    }

  cell->get_dof_indices (data.local_dof_indices);
}



template <int dim>
void
BoussinesqFlowProblem<dim>::
copy_local_to_global_stokes_system (const Assembly::CopyData::StokesSystem<dim> &data)
{
  if (rebuild_stokes_matrix == true)
    stokes_constraints.distribute_local_to_global (data.local_matrix,
						   data.local_rhs,
						   data.local_dof_indices,
						   stokes_matrix,
						   stokes_rhs);
  else
    stokes_constraints.distribute_local_to_global (data.local_rhs,
						   data.local_dof_indices,
						   stokes_rhs);
}



template <int dim>
void BoussinesqFlowProblem<dim>::assemble_stokes_system ()
{
  computing_timer.enter_section ("   Assemble Stokes system");

  if (rebuild_stokes_matrix == true)
    stokes_matrix=0;

  stokes_rhs=0;

  const QGauss<dim> quadrature_formula(stokes_degree+2);

  typedef
    FilteredIterator<typename DoFHandler<dim>::active_cell_iterator>
    SubdomainFilter;

  WorkStream::
    run (SubdomainFilter (IteratorFilters::SubdomainEqualTo
			  (Utilities::System::get_this_mpi_process(MPI_COMM_WORLD)),
			  stokes_dof_handler.begin_active()),
	 SubdomainFilter (IteratorFilters::SubdomainEqualTo
			  (Utilities::System::get_this_mpi_process(MPI_COMM_WORLD)),
			  stokes_dof_handler.end()),
	 std_cxx1x::bind (&BoussinesqFlowProblem<dim>::
			  local_assemble_stokes_system,
			  this,
			  _1,
			  _2,
			  _3),
	 std_cxx1x::bind (&BoussinesqFlowProblem<dim>::
			  copy_local_to_global_stokes_system,
			  this,
			  _1),
	 Assembly::Scratch::
	 StokesSystem<dim> (stokes_fe, quadrature_formula,
			    (update_values    |
			     update_quadrature_points  |
			     update_JxW_values |
			     (rebuild_stokes_matrix == true
			      ?
			      update_gradients
			      :
			      UpdateFlags(0))),
			    temperature_fe,
			    update_values),
	 Assembly::CopyData::
	 StokesSystem<dim> (stokes_fe));

  stokes_matrix.compress();
  stokes_rhs.compress(Add);

  rebuild_stokes_matrix = false;

  pcout << std::endl;
  computing_timer.exit_section();
}


				 // @sect5{Temperature matrix assembly}

				 // The task to be performed by the next three
				 // functions is to calculate a mass matrix
				 // and a Laplace matrix on the temperature
				 // system. These will be combined in order to
				 // yield the semi-implicit time stepping
				 // matrix that consists of the mass matrix
				 // plus a time step weight times the Laplace
				 // matrix. This function is again essentially
				 // the body of the loop over all cells from
				 // step-31.
				 //
				 // The two following functions perform
				 // similar services as the ones above.
template <int dim>
void BoussinesqFlowProblem<dim>::
local_assemble_temperature_matrix (const typename DoFHandler<dim>::active_cell_iterator &cell,
				   Assembly::Scratch::TemperatureMatrix<dim> &scratch,
				   Assembly::CopyData::TemperatureMatrix<dim> &data)
{
  const unsigned int dofs_per_cell = scratch.temperature_fe_values.get_fe().dofs_per_cell;
  const unsigned int n_q_points    = scratch.temperature_fe_values.n_quadrature_points;

  scratch.temperature_fe_values.reinit (cell);
  cell->get_dof_indices (data.local_dof_indices);

  data.local_mass_matrix = 0;
  data.local_stiffness_matrix = 0;

  for (unsigned int q=0; q<n_q_points; ++q)
    {
      for (unsigned int k=0; k<dofs_per_cell; ++k)
	{
	  scratch.grad_phi_T[k] = scratch.temperature_fe_values.shape_grad (k,q);
	  scratch.phi_T[k]      = scratch.temperature_fe_values.shape_value (k, q);
	}

      for (unsigned int i=0; i<dofs_per_cell; ++i)
	for (unsigned int j=0; j<dofs_per_cell; ++j)
	  {
	    data.local_mass_matrix(i,j)
	      += (scratch.phi_T[i] * scratch.phi_T[j]
		  *
		  scratch.temperature_fe_values.JxW(q));
	    data.local_stiffness_matrix(i,j)
	      += (EquationData::kappa * scratch.grad_phi_T[i] * scratch.grad_phi_T[j]
		  *
		  scratch.temperature_fe_values.JxW(q));
	  }
    }
}



template <int dim>
void
BoussinesqFlowProblem<dim>::
copy_local_to_global_temperature_matrix (const Assembly::CopyData::TemperatureMatrix<dim> &data)
{
  temperature_constraints.distribute_local_to_global (data.local_mass_matrix,
						      data.local_dof_indices,
						      temperature_mass_matrix);
  temperature_constraints.distribute_local_to_global (data.local_stiffness_matrix,
						      data.local_dof_indices,
						      temperature_stiffness_matrix);
}


template <int dim>
void BoussinesqFlowProblem<dim>::assemble_temperature_matrix ()
{
  if (rebuild_temperature_matrices == false)
    return;

  computing_timer.enter_section ("   Assemble temperature matrices");
  temperature_mass_matrix = 0;
  temperature_stiffness_matrix = 0;

  const QGauss<dim> quadrature_formula(temperature_degree+2);

  typedef
    FilteredIterator<typename DoFHandler<dim>::active_cell_iterator>
    SubdomainFilter;

  WorkStream::
    run (SubdomainFilter (IteratorFilters::SubdomainEqualTo
			  (Utilities::System::get_this_mpi_process(MPI_COMM_WORLD)),
			  temperature_dof_handler.begin_active()),
	 SubdomainFilter (IteratorFilters::SubdomainEqualTo
			  (Utilities::System::get_this_mpi_process(MPI_COMM_WORLD)),
			  temperature_dof_handler.end()),
	 std_cxx1x::bind (&BoussinesqFlowProblem<dim>::
			  local_assemble_temperature_matrix,
			  this,
			  _1,
			  _2,
			  _3),
	 std_cxx1x::bind (&BoussinesqFlowProblem<dim>::
			  copy_local_to_global_temperature_matrix,
			  this,
			  _1),
	 Assembly::Scratch::
	 TemperatureMatrix<dim> (temperature_fe, quadrature_formula),
	 Assembly::CopyData::
	 TemperatureMatrix<dim> (temperature_fe));

  temperature_mass_matrix.compress();
  temperature_stiffness_matrix.compress();

  rebuild_temperature_matrices = false;
  rebuild_temperature_preconditioner = true;

  computing_timer.exit_section();
}


				 // @sect5{Temperature right hand side assembly}

				 // This is the last assembly function. It
				 // calculates the right hand side of the
				 // temperature system, which includes the
				 // convection and the stabilization
				 // terms. It includes a lot of evaluations
				 // of old solutions at the quadrature
				 // points (which are necessary for
				 // calculating the artificial viscosity of
				 // stabilization), but is otherwise similar
				 // to the other assembly functions. Notice,
				 // once again, how we resolve the dilemma
				 // of having inhomogeneous boundary
				 // conditions, but just making a right hand
				 // side at this point (compare the comments
				 // for the project function): We create
				 // some matrix columns with exactly the
				 // values that would be entered for the
				 // temperature stiffness matrix, in case we
				 // have inhomogeneously constrained
				 // dofs. That will account for the correct
				 // balance of the right hand side vector
				 // with the matrix system of temperature.
template <int dim>
void BoussinesqFlowProblem<dim>::
local_assemble_temperature_rhs (const std::pair<double,double> global_T_range,
				const double                   global_max_velocity,
				const typename DoFHandler<dim>::active_cell_iterator &cell,
				Assembly::Scratch::TemperatureRHS<dim> &scratch,
				Assembly::CopyData::TemperatureRHS<dim> &data)
{
  const bool use_bdf2_scheme = (timestep_number != 0);

  const unsigned int dofs_per_cell = scratch.temperature_fe_values.get_fe().dofs_per_cell;
  const unsigned int n_q_points    = scratch.temperature_fe_values.n_quadrature_points;

  const FEValuesExtractors::Vector velocities (0);

  data.local_rhs = 0;
  data.matrix_for_bc = 0;
  cell->get_dof_indices (data.local_dof_indices);

  scratch.temperature_fe_values.reinit (cell);

  typename DoFHandler<dim>::active_cell_iterator
    stokes_cell (&triangulation,
		 cell->level(),
		 cell->index(),
		 &stokes_dof_handler);
  scratch.stokes_fe_values.reinit (stokes_cell);

  scratch.temperature_fe_values.get_function_values (old_temperature_solution,
						     scratch.old_temperature_values);
  scratch.temperature_fe_values.get_function_values (old_old_temperature_solution,
						     scratch.old_old_temperature_values);

  scratch.temperature_fe_values.get_function_gradients (old_temperature_solution,
							scratch.old_temperature_grads);
  scratch.temperature_fe_values.get_function_gradients (old_old_temperature_solution,
							scratch.old_old_temperature_grads);

  scratch.temperature_fe_values.get_function_laplacians (old_temperature_solution,
							 scratch.old_temperature_laplacians);
  scratch.temperature_fe_values.get_function_laplacians (old_old_temperature_solution,
							 scratch.old_old_temperature_laplacians);

  scratch.stokes_fe_values[velocities].get_function_values (stokes_solution,
							    scratch.old_velocity_values);
  scratch.stokes_fe_values[velocities].get_function_values (old_stokes_solution,
							    scratch.old_old_velocity_values);
  scratch.stokes_fe_values[velocities].get_function_symmetric_gradients (stokes_solution,
							       scratch.old_strain_rates);
  scratch.stokes_fe_values[velocities].get_function_symmetric_gradients (old_stokes_solution,
							       scratch.old_old_strain_rates);

  const double nu
    = compute_viscosity (scratch.old_temperature_values,
			 scratch.old_old_temperature_values,
			 scratch.old_temperature_grads,
			 scratch.old_old_temperature_grads,
			 scratch.old_temperature_laplacians,
			 scratch.old_old_temperature_laplacians,
			 scratch.old_velocity_values,
			 scratch.old_old_velocity_values,
			 scratch.old_strain_rates,
			 scratch.old_old_strain_rates,
			 global_max_velocity,
			 global_T_range.second - global_T_range.first,
			 cell->diameter());

  for (unsigned int q=0; q<n_q_points; ++q)
    {
      for (unsigned int k=0; k<dofs_per_cell; ++k)
	{
	  scratch.phi_T[k]      = scratch.temperature_fe_values.shape_value (k, q);
	  scratch.grad_phi_T[k] = scratch.temperature_fe_values.shape_grad (k,q);
	}


      const double old_Ts
	= (use_bdf2_scheme ?
	   (scratch.old_temperature_values[q] *
	    (time_step + old_time_step) / old_time_step
	    -
	    scratch.old_old_temperature_values[q] *
	    (time_step * time_step) /
	    (old_time_step * (time_step + old_time_step)))
	   :
	   scratch.old_temperature_values[q]);

      const Tensor<1,dim> ext_grad_T
	= (use_bdf2_scheme ?
	   (scratch.old_temperature_grads[q] *
	    (1+time_step/old_time_step)
	    -
	    scratch.old_old_temperature_grads[q] *
	    time_step / old_time_step)
	   :
	   scratch.old_temperature_grads[q]);

      const Tensor<1,dim> extrapolated_u
	= (use_bdf2_scheme ?
	   (scratch.old_velocity_values[q] * (1+time_step/old_time_step) -
	    scratch.old_old_velocity_values[q] * time_step/old_time_step)
	   :
	   scratch.old_velocity_values[q]);
      const SymmetricTensor<2,dim> extrapolated_strain_rate
	= (use_bdf2_scheme ?
	   (scratch.old_strain_rates[q] * (1+time_step/old_time_step) -
	    scratch.old_old_strain_rates[q] * time_step/old_time_step)
	   :
	   scratch.old_strain_rates[q]);

      const double gamma
	= ((EquationData::radiogenic_heating * EquationData::density(old_Ts) //?????? why old_Ts?
	    +
	    2 * EquationData::eta * extrapolated_strain_rate * extrapolated_strain_rate) /
	   (EquationData::density(old_Ts) * EquationData::specific_heat));
      
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	{
	  data.local_rhs(i) += (old_Ts * scratch.phi_T[i]
				-
				time_step *
				extrapolated_u * ext_grad_T * scratch.phi_T[i]
				-
				time_step *
				nu * ext_grad_T * scratch.grad_phi_T[i]
				+
				time_step *
				gamma * scratch.phi_T[i])
	                       *
	                       scratch.temperature_fe_values.JxW(q);

	  if (temperature_constraints.is_inhomogeneously_constrained(data.local_dof_indices[i]))
	    {
	      for (unsigned int j=0; j<dofs_per_cell; ++j)
		data.matrix_for_bc(j,i) += (scratch.phi_T[i] * scratch.phi_T[j] *
					    (use_bdf2_scheme ?
					     ((2*time_step + old_time_step) /
					      (time_step + old_time_step)) : 1.)
					    +
					    scratch.grad_phi_T[i] *
					    scratch.grad_phi_T[j] *
					    EquationData::kappa *
					    time_step)
		                           *
		                           scratch.temperature_fe_values.JxW(q);
	    }
	}
    }
}


template <int dim>
void
BoussinesqFlowProblem<dim>::
copy_local_to_global_temperature_rhs (const Assembly::CopyData::TemperatureRHS<dim> &data)
{
  temperature_constraints.distribute_local_to_global (data.local_rhs,
						      data.local_dof_indices,
						      temperature_rhs,
						      data.matrix_for_bc);
}



				 // In the function that runs the WorkStream
				 // for actually calculating the right hand
				 // side, we also generate the final
				 // matrix. As mentioned above, it is a sum
				 // of the mass matrix and the Laplace
				 // matrix, times some time step
				 // weight. This weight is specified by the
				 // BDF-2 time integration scheme, see the
				 // introduction in step-31. What is new in
				 // this tutorial program (in addition to
				 // the use of MPI parallelization and the
				 // WorkStream class), is that we now
				 // precompute the temperature
				 // preconditioner as well. The reason is
				 // that the setup of the IC preconditioner
				 // takes a noticable time compared to the
				 // solver because we usually only need
				 // between 10 and 20 iterations for solving
				 // the temperature system. Hence, it is
				 // more efficient to precompute the
				 // preconditioner, even though the matrix
				 // entries may slightly change because the
				 // time step might change. This is not
				 // too big a problem because we remesh every
				 // fifth time step (and regenerate the
				 // preconditioner then).
template <int dim>
void BoussinesqFlowProblem<dim>::assemble_temperature_system (const double maximal_velocity)
{
  const bool use_bdf2_scheme = (timestep_number != 0);

  if (use_bdf2_scheme == true)
    {
      temperature_matrix.copy_from (temperature_mass_matrix);
      temperature_matrix *= (2*time_step + old_time_step) /
			    (time_step + old_time_step);
      temperature_matrix.add (time_step, temperature_stiffness_matrix);
    }
  else
    {
      temperature_matrix.copy_from (temperature_mass_matrix);
      temperature_matrix.add (time_step, temperature_stiffness_matrix);
    }
  temperature_matrix.compress();

  if (rebuild_temperature_preconditioner == true)
    {
      T_preconditioner =  std_cxx1x::shared_ptr<TrilinosWrappers::PreconditionIC>
			  (new TrilinosWrappers::PreconditionIC());
      T_preconditioner->initialize (temperature_matrix);
      rebuild_temperature_preconditioner = false;
    }

  temperature_rhs = 0;

  const QGauss<dim> quadrature_formula(temperature_degree+2);
  const std::pair<double,double>
    global_T_range = get_extrapolated_temperature_range();

  typedef
    FilteredIterator<typename DoFHandler<dim>::active_cell_iterator>
    SubdomainFilter;

  WorkStream::
    run (SubdomainFilter (IteratorFilters::SubdomainEqualTo
			  (Utilities::System::get_this_mpi_process(MPI_COMM_WORLD)),
			  temperature_dof_handler.begin_active()),
	 SubdomainFilter (IteratorFilters::SubdomainEqualTo
			  (Utilities::System::get_this_mpi_process(MPI_COMM_WORLD)),
			  temperature_dof_handler.end()),
	 std_cxx1x::bind (&BoussinesqFlowProblem<dim>::
			  local_assemble_temperature_rhs,
			  this,
			  global_T_range,
			  maximal_velocity,
			  _1,
			  _2,
			  _3),
	 std_cxx1x::bind (&BoussinesqFlowProblem<dim>::
			  copy_local_to_global_temperature_rhs,
			  this,
			  _1),
	 Assembly::Scratch::
	 TemperatureRHS<dim> (temperature_fe, stokes_fe, quadrature_formula),
	 Assembly::CopyData::
	 TemperatureRHS<dim> (temperature_fe));

  temperature_rhs.compress(Add);
}




				 // @sect4{BoussinesqFlowProblem::solve}

				 // This function solves the linear systems
				 // in each time step of the Boussinesq
				 // problem. First, we
				 // work on the Stokes system and then on
				 // the temperature system. In essence, it
				 // does the same things as the respective
				 // function in step-31. However, there are
				 // a few things that we need to pay some
				 // attention to. The first thing is, as
				 // mentioned in the introduction, the way
				 // we store our solution: we keep the full
				 // vector with all degrees of freedom on
				 // each MPI node. When we enter a solver
				 // which is supposed to perform
				 // matrix-vector products with a
				 // distributed matrix, this is not the
				 // appropriate form, though. There, we will
				 // want to have the solution vector to be
				 // distributed in the same way as the
				 // matrix. So what we do first (after
				 // initializing the Schur-complement based
				 // preconditioner) is to generate a
				 // distributed vector called
				 // <code>distributed_stokes_solution</code>
				 // and put only the locally owned dofs into
				 // that, which is neatly done by the
				 // <code>operator=</code> of the Trilinos
				 // vector. Next, we need to set the
				 // pressure values at hanging nodes to
				 // zero. This we also did in step-31 in
				 // order not to disturb the Schur
				 // complement by some vector entries that
				 // actually are irrelevant during the solve
				 // stage. As a difference to step-31, here
				 // we do it only for the locally owned
				 // pressure dofs. After solving for the
				 // Stokes solution, each processor copies
				 // distributed solution back into the solution
				 // vector for which every element is locally
				 // owned.
				 //
				 // Apart from these two changes, everything
				 // is the same as in step-31, so we don't
				 // need to further comment on it.
template <int dim>
void BoussinesqFlowProblem<dim>::solve ()
{
  computing_timer.enter_section ("   Solve Stokes system");

  {
    const LinearSolvers::RightPrecond<TrilinosWrappers::PreconditionAMG,
      TrilinosWrappers::PreconditionILU>
      preconditioner (stokes_matrix, stokes_preconditioner_matrix,
		      *Mp_preconditioner, *Amg_preconditioner);

    TrilinosWrappers::MPI::BlockVector
      distributed_stokes_solution (stokes_rhs);
//    distributed_stokes_solution = stokes_solution;
    distributed_stokes_solution.block(0).reinit(stokes_solution.block(0),false,true);
    distributed_stokes_solution.block(1).reinit(stokes_solution.block(1),false,true);


    const unsigned int
      start = (distributed_stokes_solution.block(0).size() +
	       distributed_stokes_solution.block(1).local_range().first),
      end   = (distributed_stokes_solution.block(0).size() +
	       distributed_stokes_solution.block(1).local_range().second);
    for (unsigned int i=start; i<end; ++i)
      if (stokes_constraints.is_constrained (i))
	distributed_stokes_solution(i) = 0;

    SolverControl solver_control (stokes_matrix.m(), 1e-7*stokes_rhs.l2_norm());
    {
      PrimitiveVectorMemory< TrilinosWrappers::MPI::BlockVector > mem;
      SolverFGMRES<TrilinosWrappers::MPI::BlockVector>
	solver(solver_control, mem,
	       SolverFGMRES<TrilinosWrappers::MPI::BlockVector>::AdditionalData(50, true));
      solver.solve(stokes_matrix, distributed_stokes_solution, stokes_rhs,
		   preconditioner);
    }
    stokes_constraints.distribute (distributed_stokes_solution);
				     //stokes_solution = distributed_stokes_solution;
    stokes_solution.block(0).reinit(distributed_stokes_solution.block(0), false, true);
    stokes_solution.block(1).reinit(distributed_stokes_solution.block(1), false, true);

    pcout << "   "
	  << solver_control.last_step()
	  << " iterations for Stokes subsystem."
          << " reduced res by " << solver_control.last_value()/solver_control.initial_value()
	  << std::endl;
  }
  computing_timer.exit_section();


  computing_timer.enter_section ("   Assemble temperature rhs");
  {
    old_time_step = time_step;
    const double maximal_velocity = get_maximal_velocity();
    double local_time_step = 1./(1.6*dim*std::sqrt(1.*dim)) /
			     temperature_degree *
			     GridTools::minimal_cell_diameter(triangulation) /
			     std::max(1e-10,maximal_velocity);

				     // calculate the minimum allowed time step
				     // size
#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
    MPI_Allreduce (&local_time_step, &time_step, 1, MPI_DOUBLE,
		   MPI_MIN, MPI_COMM_WORLD);
#else
    time_step = local_time_step;
#endif

    pcout << "   Maximal velocity: "
	  << maximal_velocity * EquationData::year_in_seconds * 100
	  << " cm/year"
	  << std::endl;
    pcout << "   " << "Time step: "
	  << time_step/EquationData::year_in_seconds
	  << " years"
	  << std::endl;

    temperature_solution = old_temperature_solution;
    assemble_temperature_system (maximal_velocity);
  }
  computing_timer.exit_section ();

  computing_timer.enter_section ("   Solve temperature system");
  {
    SolverControl solver_control (temperature_matrix.m(),
				  1e-12*temperature_rhs.l2_norm());
    SolverCG<TrilinosWrappers::MPI::Vector>   cg (solver_control);

    TrilinosWrappers::MPI::Vector
      distributed_temperature_solution (temperature_rhs);
//    distributed_temperature_solution = temperature_solution;
    distributed_temperature_solution.reinit(temperature_solution, false, true);

    cg.solve (temperature_matrix, distributed_temperature_solution,
	      temperature_rhs, *T_preconditioner);

    temperature_constraints.distribute (distributed_temperature_solution);
//    temperature_solution = distributed_temperature_solution;
    temperature_solution.reinit(distributed_temperature_solution, false, true);

    pcout << "   "
	  << solver_control.last_step()
	  << " CG iterations for temperature" << std::endl;
    computing_timer.exit_section();

				// extract temperature range
    std::vector<double> temperature (2), global_temperature (2);
    temperature[0] = distributed_temperature_solution.trilinos_vector()[0][0],
      temperature[1] = temperature[0];
    for (unsigned int i=1; i<distributed_temperature_solution.local_size(); ++i)
      {
	temperature[0] = std::min<double> (temperature[0],
					   distributed_temperature_solution.trilinos_vector()[0][i]);
	temperature[1] = std::max<double> (temperature[1],
					   distributed_temperature_solution.trilinos_vector()[0][i]);
      }
#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
    temperature[0] *= -1.0;
    MPI_Allreduce (&temperature[0], &global_temperature[0],
		   2, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    global_temperature[0] *= -1.0;
#else
    global_temperature = local_temperature;
#endif

    pcout << "   Temperature range: "
	  << global_temperature[0] << ' ' << global_temperature[1]
	  << std::endl;
  }
}


				 // @sect4{BoussinesqFlowProblem::output_results}

template <int dim>
class BoussinesqFlowProblem<dim>::Postprocessor : public DataPostprocessor<dim>
{
  public:
    Postprocessor (const unsigned int partition);
    
    virtual
    void
    compute_derived_quantities_vector (const std::vector<Vector<double> >              &uh,
				       const std::vector<std::vector<Tensor<1,dim> > > &duh,
				       const std::vector<std::vector<Tensor<2,dim> > > &dduh,
				       const std::vector<Point<dim> >                  &normals,
				       const std::vector<Point<dim> >                  &evaluation_points,
				       std::vector<Vector<double> >                    &computed_quantities) const;

    virtual std::vector<std::string> get_names () const;

    virtual unsigned int n_output_variables() const;

    virtual
    std::vector<DataComponentInterpretation::DataComponentInterpretation>
    get_data_component_interpretation () const;

    virtual UpdateFlags get_needed_update_flags () const;

  private:
    const unsigned int partition;
};


template <int dim>
BoussinesqFlowProblem<dim>::Postprocessor::Postprocessor (const unsigned int partition)
		:
		partition (partition)
{}


template <int dim>
std::vector<std::string>
BoussinesqFlowProblem<dim>::Postprocessor::get_names() const
{
  std::vector<std::string> solution_names (dim, "velocity");
  solution_names.push_back ("p");
  solution_names.push_back ("T");
  solution_names.push_back ("friction_heating");
  solution_names.push_back ("partition");
  return solution_names;
}


template <int dim>
unsigned int
BoussinesqFlowProblem<dim>::Postprocessor::n_output_variables() const
{
  return dim + 1 + 1 + 1 + 1;
}


template <int dim>
std::vector<DataComponentInterpretation::DataComponentInterpretation>
BoussinesqFlowProblem<dim>::Postprocessor::
get_data_component_interpretation () const
{
  std::vector<DataComponentInterpretation::DataComponentInterpretation>
    interpretation (dim,
		    DataComponentInterpretation::component_is_part_of_vector);

  interpretation.push_back (DataComponentInterpretation::component_is_scalar);
  interpretation.push_back (DataComponentInterpretation::component_is_scalar);
  interpretation.push_back (DataComponentInterpretation::component_is_scalar);
  interpretation.push_back (DataComponentInterpretation::component_is_scalar);

  return interpretation;
}


template <int dim>
UpdateFlags
BoussinesqFlowProblem<dim>::Postprocessor::get_needed_update_flags() const
{
  return update_values | update_gradients;
}


template <int dim>
void
BoussinesqFlowProblem<dim>::Postprocessor::
compute_derived_quantities_vector (const std::vector<Vector<double> >              &uh,
				   const std::vector<std::vector<Tensor<1,dim> > > &duh,
				   const std::vector<std::vector<Tensor<2,dim> > > &/*dduh*/,
				   const std::vector<Point<dim> >                  &/*normals*/,
				   const std::vector<Point<dim> >                  &/*evaluation_points*/,
				   std::vector<Vector<double> >                    &computed_quantities) const
{
  const unsigned int n_quadrature_points = uh.size();
  Assert (duh.size() == n_quadrature_points,                  ExcInternalError());
  Assert (computed_quantities.size() == n_quadrature_points,  ExcInternalError());
  Assert (uh[0].size() == dim+2,                              ExcInternalError());
  Assert (computed_quantities[0].size()==n_output_variables(),ExcInternalError());

  for (unsigned int q=0; q<n_quadrature_points; ++q)
    {
				       // velocity; rescale in cm/year
      for (unsigned int d=0; d<dim; ++d)
	computed_quantities[q](d)
	  = (uh[q](d) *  EquationData::year_in_seconds * 100);

				       // pressure
      computed_quantities[q](dim) = uh[q](dim) *
				    EquationData::pressure_scaling;
  
				       // temperature
      computed_quantities[q](dim+1) = uh[q](dim+1);

				       // friction heating
      Tensor<2,dim> grad_u;
      for (unsigned int d=0; d<dim; ++d)
	grad_u[d] = duh[q][d];
      const SymmetricTensor<2,dim> strain_rate = symmetrize (grad_u);
      computed_quantities[q](dim+2) = 2 * EquationData::eta * strain_rate * strain_rate;

      computed_quantities[q](dim+3) = partition;
    }
}


				 // This function does mostly what the
				 // corresponding one did in to
				 // step-31, in particular merging
				 // data from the two DoFHandler
				 // objects (for the Stokes and the
				 // temperature parts of the problem)
				 // into one is the same. There are
				 // three minor changes: we make sure
				 // that only a single processor
				 // actually does some work here; take
				 // care of scaling variables in a
				 // useful way; and in addition to the
				 // Stokes and temperature parts in
				 // the <code>joint_fe</code> finite
				 // element, we also add a piecewise
				 // constant field that denotes the
				 // subdomain id a cell corresponds
				 // to. This allows us to visualize
				 // the partitioning of the domain. As
				 // a consequence, we also have to
				 // change the assertion about the
				 // number of degrees of freedom in
				 // the joint DoFHandler object (which
				 // is now equal to the number of
				 // Stokes degrees of freedom plus the
				 // temperature degrees of freedom
				 // plus the number of active cells as
				 // that is the number of partition
				 // variables we want to add), and
				 // adjust the number of elements in
				 // the arrays we use to name the
				 // components of the joint solution
				 // vector and to identify which of
				 // these components are scalars or
				 // parts of dim-dimensional vectors.
				 //
				 // As for scaling: as mentioned in
				 // the introduction, to keep the
				 // Stokes equations properly scaled
				 // and symmetric, we introduced a new
				 // pressure $\hat p =
				 // \frac{L}{\eta}p$. What we really
				 // wanted, however, was the original
				 // pressure $p$, so while copying
				 // data from the Stokes DoFHandler
				 // into the joint one, we undo this
				 // scaling. While we're at it messing
				 // with the results of the
				 // simulation, we do two more things:
				 // First, the pressure is only
				 // defined up to a constant. To make
				 // it more easily comparable, we
				 // compute the minimal value of the
				 // pressure computed and shift all
				 // values up by that amount -- in
				 // essence making all pressure
				 // variables positive or
				 // zero. Secondly, let's also take
				 // care of the awkward units we use
				 // for the velocity: it is computed
				 // in SI units of meters per second,
				 // which of course is a very small
				 // number in the earth mantle. We
				 // therefore rescale things into
				 // centimeters per year, the unit
				 // commonly used in geophysics.
template <int dim>
void BoussinesqFlowProblem<dim>::output_results ()
{
  computing_timer.enter_section ("Postprocessing");
  Vector<float> estimated_error_per_cell (triangulation.n_active_cells());

  KellyErrorEstimator<dim>::estimate (temperature_dof_handler,
                                      QGauss<dim-1>(temperature_degree+1),
                                      typename FunctionMap<dim>::type(),
                                      temperature_solution,
                                      estimated_error_per_cell);

  const FESystem<dim> joint_fe (stokes_fe, 1,
                                temperature_fe, 1);

  DoFHandler<dim> joint_dof_handler (triangulation);
  joint_dof_handler.distribute_dofs (joint_fe);
  Assert (joint_dof_handler.n_dofs() ==
	  stokes_dof_handler.n_dofs() + temperature_dof_handler.n_dofs(),
	  ExcInternalError());

  TrilinosWrappers::MPI::Vector joint_solution;
  joint_solution.reinit (joint_dof_handler.locally_owned_dofs(), MPI_COMM_WORLD);
  
  {
    //double minimal_pressure = stokes_solution.block(1)(0);
    //for (unsigned int i=0; i<stokes_solution.block(1).size(); ++i)
    //  minimal_pressure = std::min<double> (stokes_solution.block(1)(i),
    //					   minimal_pressure);

    std::vector<unsigned int> local_joint_dof_indices (joint_fe.dofs_per_cell);
    std::vector<unsigned int> local_stokes_dof_indices (stokes_fe.dofs_per_cell);
    std::vector<unsigned int> local_temperature_dof_indices (temperature_fe.dofs_per_cell);

    typename DoFHandler<dim>::active_cell_iterator
      joint_cell       = joint_dof_handler.begin_active(),
      joint_endc       = joint_dof_handler.end(),
      stokes_cell      = stokes_dof_handler.begin_active(),
      temperature_cell = temperature_dof_handler.begin_active();
    for (; joint_cell!=joint_endc;
	 ++joint_cell, ++stokes_cell, ++temperature_cell)
      if (!joint_cell->is_artificial() && !joint_cell->is_ghost())
	{
	  joint_cell->get_dof_indices (local_joint_dof_indices);
	  stokes_cell->get_dof_indices (local_stokes_dof_indices);
	  temperature_cell->get_dof_indices (local_temperature_dof_indices);

	  for (unsigned int i=0; i<joint_fe.dofs_per_cell; ++i)
	    if (joint_fe.system_to_base_index(i).first.first == 0)
	      {
		Assert (joint_fe.system_to_base_index(i).second
			<
			local_stokes_dof_indices.size(),
			ExcInternalError());

		joint_solution(local_joint_dof_indices[i])
		  = stokes_solution(local_stokes_dof_indices
				    [joint_fe.system_to_base_index(i).second]);
	      }
	    else
	      {
		Assert (joint_fe.system_to_base_index(i).first.first == 1,
			ExcInternalError());
		Assert (joint_fe.system_to_base_index(i).second
			<
			local_temperature_dof_indices.size(),
			ExcInternalError());
		joint_solution(local_joint_dof_indices[i])
		  = temperature_solution(local_temperature_dof_indices
					 [joint_fe.system_to_base_index(i).second]);
	      }
	}
  }


  IndexSet locally_relevant_joint_dofs(joint_dof_handler.n_dofs());
  DoFTools::extract_locally_relevant_dofs (joint_dof_handler, locally_relevant_joint_dofs);
  TrilinosWrappers::MPI::Vector locally_relevant_joint_solution;
  locally_relevant_joint_solution.reinit (locally_relevant_joint_dofs, MPI_COMM_WORLD);
  locally_relevant_joint_solution = joint_solution;

  Postprocessor postprocessor (Utilities::System::
			       get_this_mpi_process(MPI_COMM_WORLD));

  DataOut<dim> data_out;
  data_out.attach_dof_handler (joint_dof_handler);
  data_out.add_data_vector (locally_relevant_joint_solution, postprocessor);
  data_out.build_patches ();

  static int out_index=0;
  const std::string filename = ("solution-" +
				Utilities::int_to_string (out_index, 5) +
				"." +
				Utilities::int_to_string
				(triangulation.locally_owned_subdomain(), 4) +
				".vtu");
  std::ofstream output (filename.c_str());
  data_out.write_vtu (output);

  if (Utilities::System::get_this_mpi_process(MPI_COMM_WORLD) == 0)
    {
      std::vector<std::string> filenames;
      for (unsigned int i=0; i<Utilities::System::get_n_mpi_processes(MPI_COMM_WORLD); ++i)
	filenames.push_back (std::string("solution-") +
			     Utilities::int_to_string (out_index, 5) +
			     "." +
			     Utilities::int_to_string(i, 4) +
			     ".vtu");
      const std::string
	master_filename = ("solution-" +
			   Utilities::int_to_string (out_index, 5) +
			   ".pvtu");
      std::ofstream master (master_filename.c_str());
      data_out.write_pvtu_record (master, filenames);
    }

  computing_timer.exit_section ();
  out_index++;
}



				 // @sect4{BoussinesqFlowProblem::refine_mesh}

				 // This function isn't really new
				 // either. Since the
				 // <code>setup_dofs</code> function
				 // that we call in the middle has its
				 // own timer section, we split timing
				 // this function into two
				 // sections. It will also allow us to
				 // easily identify which of the two
				 // is more expensive.
				 //
				 // One thing of note, however, is
				 // that we don't want to compute all
				 // error indicators on all cells, of
				 // course. Rather, it would be nice
				 // if each processor could only
				 // compute the error indicators for
				 // those cells it actually
				 // owns. However, in order for mesh
				 // refinement to proceed in the same
				 // way on all processors, all
				 // processors would have to exchange
				 // their refinement indicators. We do
				 // so in two steps: first, we call
				 // the KellyErrorEstimator::estimate
				 // function with an argument (usually
				 // defaulted, but explicitly given
				 // here) thatindicates the subdomain
				 // id of all those cells that we want
				 // to work on; note that this means
				 // that we also have to specify
				 // values for all those default
				 // arguments that lie before the one
				 // we want to give.
				 //
				 // Secondly, we need to exchange the
				 // data. To do this, we could add up
				 // the refinement indicators from all
				 // processors, since they all only
				 // worked on a disjoint subset of the
				 // elements of the vector that holds
				 // these indicators. We could set up
				 // a distributed Trilinos vector for
				 // this, but that appears
				 // unnecessarily complicated because
				 // we would have to specify a
				 // partition of this vector, and none
				 // appears immediately
				 // obvious. Rather, we want to use
				 // the Trilinos communicator class to
				 // this for us, taking the local
				 // indicators as a collection of
				 // floating point values rather than
				 // a linear algebra
				 // vector. Unfortunately, the
				 // Trilinos communicator class
				 // doesn't appear to have function
				 // that wraps around the MPI add
				 // function; it has one that computes
				 // the maximum of a bunch of values,
				 // though, which in our case is
				 // equally good -- maybe even better,
				 // in case two processors should
				 // compute values for the same cell
				 // (which they shouldn't of course,
				 // unless we have made a mistake in
				 // specifying the arguments to the
				 // estimate function below). There is
				 // little snag again, however, that
				 // makes this a bit awkward: the
				 // Trilinos communicator class can
				 // take the maximum over all
				 // processors for each element of a
				 // vector, but only if the vector
				 // contains doubles. The vector
				 // returned by the
				 // KellyErrorEstimator::estimate
				 // function, on the other hand, has
				 // floats as its data type. An ugly,
				 // if workable way, is therefore to
				 // compute the indicators as floats,
				 // convert the vector to doubles, and
				 // form the maximum of that.
				 //
				 // At the end of this chain of
				 // events, every processors has the
				 // complete set of refinement
				 // indicators, and the rest of the
				 // function proceeds as before.
template <int dim>
void BoussinesqFlowProblem<dim>::refine_mesh (const unsigned int max_grid_level)
{
  computing_timer.enter_section ("Refine mesh structure, part 1");
  Vector<float> estimated_error_per_cell (triangulation.n_active_cells());

  KellyErrorEstimator<dim>::estimate (temperature_dof_handler,
				      QGauss<dim-1>(temperature_degree+1),
				      typename FunctionMap<dim>::type(),
				      temperature_solution,
 				      estimated_error_per_cell,
				      std::vector<bool>(),
				      0,
				      0,
				      triangulation.locally_owned_subdomain());

  parallel::distributed::GridRefinement::
    refine_and_coarsen_fixed_fraction (triangulation,
				       estimated_error_per_cell,
				       0.6, 0.2);

				   // limit maximum refinement level
  if (triangulation.n_levels() > max_grid_level)
    for (typename Triangulation<dim>::active_cell_iterator
	   cell = triangulation.begin_active(max_grid_level);
	 cell != triangulation.end(); ++cell)
      cell->clear_refine_flag ();

  std::vector<const TrilinosWrappers::MPI::Vector*> x_temperature (2);
  x_temperature[0] = &temperature_solution;
  x_temperature[1] = &old_temperature_solution;
  std::vector<const TrilinosWrappers::MPI::BlockVector*> x_stokes (2);
  x_stokes[0] = &stokes_solution;
  x_stokes[1] = &old_stokes_solution;

  parallel::distributed::SolutionTransfer<dim,TrilinosWrappers::MPI::Vector>
    temperature_trans(temperature_dof_handler);
  parallel::distributed::SolutionTransfer<dim,TrilinosWrappers::MPI::BlockVector>
    stokes_trans(stokes_dof_handler);

  triangulation.prepare_coarsening_and_refinement();
  temperature_trans.prepare_for_coarsening_and_refinement(x_temperature);
  stokes_trans.prepare_for_coarsening_and_refinement(x_stokes);

  triangulation.execute_coarsening_and_refinement ();
  computing_timer.exit_section();

  setup_dofs ();

  computing_timer.enter_section ("Refine mesh structure, part 2");

  {
    TrilinosWrappers::MPI::Vector
      distributed_temp1 (temperature_rhs);
    TrilinosWrappers::MPI::Vector
      distributed_temp2 (temperature_rhs);

    std::vector<TrilinosWrappers::MPI::Vector*> tmp (2);
    tmp[0] = &(distributed_temp1);
    tmp[1] = &(distributed_temp2);
    temperature_trans.interpolate(tmp);

    //  temperature_solution = distributed_temp1;
    temperature_solution.reinit(distributed_temp1, false, true);
    //  old_temperature_solution = distributed_temp2;
    old_temperature_solution.reinit(distributed_temp2, false, true);
  }

  {
    TrilinosWrappers::MPI::BlockVector
      distributed_stokes (stokes_rhs);
    TrilinosWrappers::MPI::BlockVector
      old_distributed_stokes (stokes_rhs);
    std::vector<TrilinosWrappers::MPI::BlockVector*> stokes_tmp (2);
    stokes_tmp[0] = &(distributed_stokes);
    stokes_tmp[1] = &(old_distributed_stokes);

    stokes_trans.interpolate (stokes_tmp);
    //  stokes_solution = distributed_stokes;
    stokes_solution.block(0).reinit(distributed_stokes.block(0), false, true);
    stokes_solution.block(1).reinit(distributed_stokes.block(1), false, true);
    old_stokes_solution.block(0).reinit(old_distributed_stokes.block(0), false, true);
    old_stokes_solution.block(1).reinit(old_distributed_stokes.block(1), false, true);
  }

  computing_timer.exit_section();
}



				 // @sect4{BoussinesqFlowProblem::run}

				 // This is the final function in this
				 // class. It actually runs the program. It
				 // is, once more, very similar to
				 // step-31. The only thing that really
				 // changed is that we use the
				 // <code>project_temperature_field()</code>
				 // function instead of the library function
				 // <code>VectorTools::project</code>, the
				 // rest is as before.
template <int dim>
void BoussinesqFlowProblem<dim>::run (unsigned int ref)
{
  pcout << "this is step-32. ref=" << ref << std::endl;
  const unsigned int initial_refinement = ref;//(dim == 2 ? 5 : 2);
  const unsigned int n_pre_refinement_steps = 2;//(dim == 2 ? 4 : 2);

  GridGenerator::hyper_shell (triangulation,
			      Point<dim>(),
			      EquationData::R0,
			      EquationData::R1,
			      12,
			      true);
  static HyperShellBoundary<dim> boundary;
  triangulation.set_boundary (0, boundary);
  triangulation.set_boundary (1, boundary);

				   //GridGenerator::hyper_cube (triangulation, EquationData::R0, EquationData::R1);

  global_Omega_diameter = GridTools::diameter (triangulation);

  triangulation.refine_global (initial_refinement);

  setup_dofs();

  unsigned int pre_refinement_step = 0;

  start_time_iteration:

  project_temperature_field ();

  timestep_number           = 0;
  time_step = old_time_step = 0;

  double time = 0;

  do
    {
      pcout << "Timestep " << timestep_number
	    << ":  t=" << time/EquationData::year_in_seconds
	    << " years"
	    << std::endl;

      assemble_stokes_system ();
      build_stokes_preconditioner ();
      assemble_temperature_matrix ();

      solve ();

      pcout << std::endl;

      if ((timestep_number == 0) &&
	  (pre_refinement_step < n_pre_refinement_steps))
	{
	  refine_mesh (initial_refinement + n_pre_refinement_steps);
	  ++pre_refinement_step;
	  goto start_time_iteration;
	}
      else
	if ((timestep_number > 0) && (timestep_number % 10 == 0))
	  refine_mesh (initial_refinement + n_pre_refinement_steps);

      if (timestep_number % 50 == 0 &&
	  Utilities::System::get_n_mpi_processes(MPI_COMM_WORLD) <= 100)
	output_results ();

      time += time_step;
      ++timestep_number;

      TrilinosWrappers::MPI::BlockVector old_old_stokes_solution;
      old_old_stokes_solution      = old_stokes_solution;
      old_stokes_solution          = stokes_solution;
      old_old_temperature_solution = old_temperature_solution;
      old_temperature_solution     = temperature_solution;
      if (old_time_step > 0)
      	{
      	  stokes_solution.sadd (1.+time_step/old_time_step, -time_step/old_time_step,
      				old_old_stokes_solution);
      	  temperature_solution.sadd (1.+time_step/old_time_step, 
      				     -time_step/old_time_step,
      				     old_old_temperature_solution);
      	}
    }
  while (time <= EquationData::end_time);
}



				 // @sect3{The <code>main</code> function}

				 // This is copied verbatim from step-31:
int main (int argc, char *argv[])
{
  try
    {
      deallog.depth_console (0);

      Utilities::System::MPI_InitFinalize mpi_initialization(argc, argv);

      const int dim = 3;
      BoussinesqFlowProblem<dim> flow_problem;

      unsigned int ref = (dim == 2 ? 5 : 2);
      if (argc>=2)
	ref = (unsigned int)Utilities::string_to_int(argv[1]);

      flow_problem.run (ref);
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;

      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }

  return 0;
}
