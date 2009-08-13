/* $Id$ */
/* Author: Martin Kronbichler, Uppsala University,
           Wolfgang Bangerth, Texas A&M University 2007, 2008, 2009 */
/*                                                                */
/*    Copyright (C) 2008, 2009 by the deal.II authors */
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
#include <lac/constraint_matrix.h>
#include <lac/block_sparsity_pattern.h>
#include <lac/trilinos_block_vector.h>
#include <lac/trilinos_sparse_matrix.h>
#include <lac/trilinos_block_sparse_matrix.h>
#include <lac/trilinos_precondition.h>

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
#include <sstream>

				 // This is the only include file that is
				 // new: We use Trilinos for defining the
				 // %parallel partitioning of the matrices
				 // and vectors, and as explained in the
				 // introduction, an <code>Epetra_Map</code>
				 // is the Trilinos data structure for the
				 // definition of which part of a
				 // distributed vector is stored locally:
#include <Epetra_Map.h>


				 // Next, we import all deal.II
				 // names into global namespace:
using namespace dealii;

				 // @sect3{Equation data}

				 // This program is mainly an extension of
				 // step-31 to operate in %parallel, so the
				 // equation data remains the same.
namespace EquationData
{
  const double eta = 1;
  const double kappa = 1e-6;
  const double Rayleigh_number = 10;



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
  TemperatureInitialValues<dim>::value (const Point<dim>  &,
					const unsigned int) const
  {
				     /* Data for shell problem */
				     /*return (p.norm() < 0.55+0.02*std::sin(p[0]*20) ? 1 : 0);*/

				     /* Data for cube problem */
    return 0.;
  }


  template <int dim>
  void
  TemperatureInitialValues<dim>::vector_value (const Point<dim> &p,
					       Vector<double>   &values) const
  {
    for (unsigned int c=0; c<this->n_components; ++c)
      values(c) = TemperatureInitialValues<dim>::value (p, c);
  }



  template <int dim>
  class TemperatureRightHandSide : public Function<dim>
  {
    public:
      TemperatureRightHandSide () : Function<dim>(1) {}

      virtual double value (const Point<dim>   &p,
			    const unsigned int  component = 0) const;

      virtual void vector_value (const Point<dim> &p,
				 Vector<double>   &value) const;
  };



  template <int dim>
  double
  TemperatureRightHandSide<dim>::value (const Point<dim>  &p,
					const unsigned int component) const
  {
				     /* Data for shell problem. */
				     /*    return 0; */

				     /* Data for cube problem. */
    Assert (component == 0,
	    ExcMessage ("Invalid operation for a scalar function."));

    Assert ((dim==2) || (dim==3), ExcNotImplemented());

    static const Point<dim> source_centers[3]
      = { (dim == 2 ? Point<dim>(.3,.1) : Point<dim>(.3,.5,.1)),
	  (dim == 2 ? Point<dim>(.45,.1) : Point<dim>(.45,.5,.1)),
	  (dim == 2 ? Point<dim>(.75,.1) : Point<dim>(.75,.5,.1)) };
    static const double source_radius
      = (dim == 2 ? 1./32 : 1./8);

    return ((source_centers[0].distance (p) < source_radius)
	    ||
	    (source_centers[1].distance (p) < source_radius)
	    ||
	    (source_centers[2].distance (p) < source_radius)
	    ?
	    1
	    :
	    0);
  }


  template <int dim>
  void
  TemperatureRightHandSide<dim>::vector_value (const Point<dim> &p,
					       Vector<double>   &values) const
  {
    for (unsigned int c=0; c<this->n_components; ++c)
      values(c) = TemperatureRightHandSide<dim>::value (p, c);
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
  class BlockSchurPreconditioner : public Subscriptor
  {
    public:
      BlockSchurPreconditioner (
	const TrilinosWrappers::BlockSparseMatrix  &S,
	const PreconditionerMp                     &Mppreconditioner,
	const PreconditionerA                      &Apreconditioner);

      void vmult (TrilinosWrappers::MPI::BlockVector       &dst,
		  const TrilinosWrappers::MPI::BlockVector &src) const;

    private:
      const SmartPointer<const TrilinosWrappers::BlockSparseMatrix> stokes_matrix;
      const PreconditionerMp &mp_preconditioner;
      const PreconditionerA  &a_preconditioner;
      mutable TrilinosWrappers::MPI::Vector tmp;
  };



  template <class PreconditionerA, class PreconditionerMp>
  BlockSchurPreconditioner<PreconditionerA, PreconditionerMp>::
  BlockSchurPreconditioner(const TrilinosWrappers::BlockSparseMatrix &S,
			   const PreconditionerMp                    &Mppreconditioner,
			   const PreconditionerA                     &Apreconditioner)
		  :
		  stokes_matrix     (&S),
		  mp_preconditioner (Mppreconditioner),
		  a_preconditioner  (Apreconditioner),
		  tmp               (stokes_matrix->block(1,1).row_partitioner())
  {}



  template <class PreconditionerA, class PreconditionerMp>
  void BlockSchurPreconditioner<PreconditionerA, PreconditionerMp>::vmult (
    TrilinosWrappers::MPI::BlockVector       &dst,
    const TrilinosWrappers::MPI::BlockVector &src) const
  {
    a_preconditioner.vmult (dst.block(0), src.block(0));
    stokes_matrix->block(1,0).residual(tmp, dst.block(0), src.block(1));
    tmp *= -1;
    mp_preconditioner.vmult (dst.block(1), tmp);
  }
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
		    old_temperature_values (stokes_quadrature.n_quadrature_points)
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

	std::vector<double>         old_temperature_values;
	std::vector<double>         old_old_temperature_values;
	std::vector<Tensor<1,dim> > old_temperature_grads;
	std::vector<Tensor<1,dim> > old_old_temperature_grads;
	std::vector<double>         old_temperature_laplacians;
	std::vector<double>         old_old_temperature_laplacians;

	std::vector<double>         gamma_values;
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
				      update_values),
		    phi_T (temperature_fe.dofs_per_cell),
		    grad_phi_T (temperature_fe.dofs_per_cell),

		    old_velocity_values (quadrature.n_quadrature_points),
		    old_old_velocity_values (quadrature.n_quadrature_points),

		    old_temperature_values (quadrature.n_quadrature_points),
		    old_old_temperature_values(quadrature.n_quadrature_points),
		    old_temperature_grads(quadrature.n_quadrature_points),
		    old_old_temperature_grads(quadrature.n_quadrature_points),
		    old_temperature_laplacians(quadrature.n_quadrature_points),
		    old_old_temperature_laplacians(quadrature.n_quadrature_points),

		    gamma_values (quadrature.n_quadrature_points)
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

		    old_temperature_values (scratch.old_temperature_values),
		    old_old_temperature_values (scratch.old_old_temperature_values),
		    old_temperature_grads (scratch.old_temperature_grads),
		    old_old_temperature_grads (scratch.old_old_temperature_grads),
		    old_temperature_laplacians (scratch.old_temperature_laplacians),
		    old_old_temperature_laplacians (scratch.old_old_temperature_laplacians),

		    gamma_values (scratch.gamma_values)
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
    };

    template <int dim>
    TemperatureRHS<dim>::
    TemperatureRHS (const FiniteElement<dim> &temperature_fe)
		    :
		    local_rhs (temperature_fe.dofs_per_cell),
		    local_dof_indices (temperature_fe.dofs_per_cell)
    {}


    template <int dim>
    TemperatureRHS<dim>::
    TemperatureRHS (const TemperatureRHS &data)
		    :
		    local_rhs (data.local_rhs),
		    local_dof_indices (data.local_dof_indices)
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
				 // Moreover, we include an MPI communicator
				 // and an Epetra_Map (see the introduction)
				 // that are needed for communication and
				 // data exchange if the Trilinos matrices
				 // and vectors are distributed over several
				 // processors. Finally, the
				 // <code>pcout</code> (for <i>%parallel
				 // <code>std::cout</code></i>) object is
				 // used to simplify writing output: each
				 // MPI process can use this to generate
				 // output as usual, but since each of these
				 // processes will produce the same output
				 // it will just be replicated many times
				 // over; with the ConditionalOStream class,
				 // only the output generated by one task
				 // will actually be printed to screen,
				 // whereas the output by all the other
				 // threads will simply be forgotten.
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
				 // actually be considered parallel
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
    void run ();

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
		      const std::vector<double>          &gamma_values,
		      const double                        global_u_infty,
		      const double                        global_T_variation,
		      const double                        cell_diameter) const;


    const Epetra_Comm                  &trilinos_communicator;

    ConditionalOStream                  pcout;

    Triangulation<dim>                  triangulation;
    double                              global_Omega_diameter;

    const unsigned int                  stokes_degree;
    FESystem<dim>                       stokes_fe;
    DoFHandler<dim>                     stokes_dof_handler;
    ConstraintMatrix                    stokes_constraints;

    std::vector<Epetra_Map>             stokes_partitioner;
    TrilinosWrappers::BlockSparseMatrix stokes_matrix;
    TrilinosWrappers::BlockSparseMatrix stokes_preconditioner_matrix;

    TrilinosWrappers::BlockVector       stokes_solution;
    TrilinosWrappers::BlockVector       old_stokes_solution;
    TrilinosWrappers::MPI::BlockVector  stokes_rhs;


    const unsigned int                  temperature_degree;
    FE_Q<dim>                           temperature_fe;
    DoFHandler<dim>                     temperature_dof_handler;
    ConstraintMatrix                    temperature_constraints;

    Epetra_Map                          temperature_partitioner;
    TrilinosWrappers::SparseMatrix      temperature_mass_matrix;
    TrilinosWrappers::SparseMatrix      temperature_stiffness_matrix;
    TrilinosWrappers::SparseMatrix      temperature_matrix;

    TrilinosWrappers::Vector            temperature_solution;
    TrilinosWrappers::Vector            old_temperature_solution;
    TrilinosWrappers::Vector            old_old_temperature_solution;
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

    void setup_stokes_matrix ();
    void setup_stokes_preconditioner ();
    void setup_temperature_matrices ();

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
                trilinos_communicator (Utilities::Trilinos::comm_world()),
		pcout (std::cout,
		       (Utilities::Trilinos::
			get_this_mpi_process(trilinos_communicator)
			== 0)),

		triangulation (Triangulation<dim>::maximum_smoothing),

                stokes_degree (1),
                stokes_fe (FE_Q<dim>(stokes_degree+1), dim,
			   FE_Q<dim>(stokes_degree), 1),
		stokes_dof_handler (triangulation),

		temperature_degree (2),
		temperature_fe (temperature_degree),
                temperature_dof_handler (triangulation),

		temperature_partitioner (0, 0, trilinos_communicator),

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
				 // function of the Trilinos
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
	Utilities::Trilinos::get_this_mpi_process(trilinos_communicator))
      {
	fe_values.reinit (cell);
	fe_values[velocities].get_function_values (stokes_solution,
						   velocity_values);

	for (unsigned int q=0; q<n_q_points; ++q)
	  max_local_velocity = std::max (max_local_velocity,
					 velocity_values[q].norm());
      }

  double max_velocity = 0.;
  trilinos_communicator.MaxAll(&max_local_velocity, &max_velocity, 1);

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

  if (timestep_number != 0)
    {
      double min_local_temperature = (1. + time_step/old_time_step) *
	                             old_temperature_solution.linfty_norm()
	                             +
			             time_step/old_time_step *
			             old_old_temperature_solution.linfty_norm(),
	     max_local_temperature = -min_local_temperature;

      typename DoFHandler<dim>::active_cell_iterator
	cell = temperature_dof_handler.begin_active(),
	endc = temperature_dof_handler.end();
      for (; cell!=endc; ++cell)
	if (cell->subdomain_id() ==
	    Utilities::Trilinos::get_this_mpi_process(trilinos_communicator))
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

      double min_temperature, max_temperature;

      trilinos_communicator.MaxAll(&max_local_temperature, &max_temperature, 1);
      trilinos_communicator.MinAll(&min_local_temperature, &min_temperature, 1);

      return std::make_pair(min_temperature, max_temperature);
    }
  else
    {
      double min_local_temperature = old_temperature_solution.linfty_norm(),
	     max_local_temperature = -min_local_temperature;

      typename DoFHandler<dim>::active_cell_iterator
	cell = temperature_dof_handler.begin_active(),
	endc = temperature_dof_handler.end();
      for (; cell!=endc; ++cell)
	if (cell->subdomain_id() ==
	    Utilities::Trilinos::get_this_mpi_process(trilinos_communicator))
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

      double min_temperature, max_temperature;

      trilinos_communicator.MaxAll(&max_local_temperature, &max_temperature, 1);
      trilinos_communicator.MinAll(&min_local_temperature, &min_temperature, 1);

      return std::make_pair(min_temperature, max_temperature);
    }
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
		   const std::vector<double>          &gamma_values,
		   const double                        global_u_infty,
		   const double                        global_T_variation,
		   const double                        cell_diameter) const
{
  const double beta = 0.015 * dim;
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

      const double dT_dt = (old_temperature[q] - old_old_temperature[q])
			   / old_time_step;
      const double u_grad_T = u * (old_temperature_grads[q] +
				   old_old_temperature_grads[q]) / 2;

      const double kappa_Delta_T = EquationData::kappa
				   * (old_temperature_laplacians[q] +
				      old_old_temperature_laplacians[q]) / 2;

      const double residual
	= std::abs((dT_dt + u_grad_T - kappa_Delta_T - gamma_values[q]) *
		   std::pow((old_temperature[q]+old_old_temperature[q]) / 2,
			    alpha-1.));

      max_residual = std::max (residual,        max_residual);
      max_velocity = std::max (std::sqrt (u*u), max_velocity);
    }

  const double global_scaling = global_u_infty * global_T_variation /
				std::pow(global_Omega_diameter, alpha - 2.);

  return (beta *
	  max_velocity *
	  std::min (cell_diameter,
		    std::pow(cell_diameter,alpha) *
		    max_residual / global_scaling));
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

  std::vector<unsigned int> dofs (dofs_per_cell);
  Vector<double> cell_vector (dofs_per_cell);

  typename DoFHandler<dim>::active_cell_iterator
    cell = temperature_dof_handler.begin_active(),
    endc = temperature_dof_handler.end();

  std::vector<double> rhs_values(n_q_points);

  TrilinosWrappers::MPI::Vector
    rhs (temperature_mass_matrix.row_partitioner()),
    solution (temperature_mass_matrix.row_partitioner());

  for (; cell!=endc; ++cell)
    if (cell->subdomain_id() ==
	Utilities::Trilinos::get_this_mpi_process(trilinos_communicator))
      {
	fe_values.reinit(cell);

	const std::vector<double> &weights   = fe_values.get_JxW_values ();
	EquationData::TemperatureInitialValues<dim>().value_list
	  (fe_values.get_quadrature_points(), rhs_values);

	cell_vector = 0;
	for (unsigned int point=0; point<n_q_points; ++point)
	  for (unsigned int i=0; i<dofs_per_cell; ++i)
	    cell_vector(i) += rhs_values[point] *
			      fe_values.shape_value(i,point) *
			      weights[point];

	cell->get_dof_indices (dofs);

	temperature_constraints.distribute_local_to_global (cell_vector,
							    dofs,
							    rhs);
      }

  ReductionControl  control(5*rhs.size(), 0., 1e-12, false, false);
  GrowingVectorMemory<TrilinosWrappers::MPI::Vector> memory;
  SolverCG<TrilinosWrappers::MPI::Vector> cg(control,memory);

  TrilinosWrappers::PreconditionIC preconditioner_mass;
  preconditioner_mass.initialize(temperature_mass_matrix);

  cg.solve (temperature_mass_matrix, solution, rhs, preconditioner_mass);

  old_temperature_solution = solution;
  temperature_constraints.distribute (old_temperature_solution);
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
				 // can easily be run in parallel on
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
void BoussinesqFlowProblem<dim>::setup_stokes_matrix ()
{
  stokes_matrix.clear ();

  TrilinosWrappers::BlockSparsityPattern sp (stokes_partitioner);

  Table<2,DoFTools::Coupling> coupling (dim+1, dim+1);

  for (unsigned int c=0; c<dim+1; ++c)
    for (unsigned int d=0; d<dim+1; ++d)
      if (! ((c==dim) && (d==dim)))
	coupling[c][d] = DoFTools::always;
      else
	coupling[c][d] = DoFTools::none;

  DoFTools::make_sparsity_pattern (stokes_dof_handler, coupling, sp,
				   stokes_constraints, false,
				   Utilities::Trilinos::
				   get_this_mpi_process(trilinos_communicator));
  sp.compress();

  stokes_matrix.reinit (sp);
}



template <int dim>
void BoussinesqFlowProblem<dim>::setup_stokes_preconditioner ()
{
  Amg_preconditioner.reset ();
  Mp_preconditioner.reset ();

  stokes_preconditioner_matrix.clear ();

  TrilinosWrappers::BlockSparsityPattern sp (stokes_partitioner);

  Table<2,DoFTools::Coupling> coupling (dim+1, dim+1);
  for (unsigned int c=0; c<dim+1; ++c)
    for (unsigned int d=0; d<dim+1; ++d)
      if (c == d)
	coupling[c][d] = DoFTools::always;
      else
	coupling[c][d] = DoFTools::none;

  DoFTools::make_sparsity_pattern (stokes_dof_handler, coupling, sp,
				   stokes_constraints, false,
				   Utilities::Trilinos::
				   get_this_mpi_process(trilinos_communicator));
  sp.compress();

  stokes_preconditioner_matrix.reinit (sp);
}


template <int dim>
void BoussinesqFlowProblem<dim>::setup_temperature_matrices ()
{
  T_preconditioner.reset ();
  temperature_mass_matrix.clear ();
  temperature_stiffness_matrix.clear ();
  temperature_matrix.clear ();

  TrilinosWrappers::SparsityPattern sp (temperature_partitioner);
  DoFTools::make_sparsity_pattern (temperature_dof_handler, sp,
				   temperature_constraints, false,
				   Utilities::Trilinos::
				   get_this_mpi_process(trilinos_communicator));
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
				 // work on.
				 //
				 // After this, we have to set up the
				 // various partitioners (of type
				 // <code>Epetra_Map</code>, see the
				 // introduction) that describe which parts
				 // of each matrix or vector will be stored
				 // where, then call the functions that
				 // actually set up the matrices
				 // (concurrently if on a single processor,
				 // but sequentially if we need MPI
				 // communications), and at the end also
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
				 // (Utilities:: Trilinos::
				 // get_this_mpi_process
				 // (trilinos_communicator) == 0)</code>,
				 // hardly a pretty solution.
template <int dim>
void BoussinesqFlowProblem<dim>::setup_dofs ()
{
  computing_timer.enter_section("Setup dof systems");
  std::vector<unsigned int> stokes_sub_blocks (dim+1,0);
  stokes_sub_blocks[dim] = 1;

  GridTools::partition_triangulation (Utilities::Trilinos::
				      get_n_mpi_processes(trilinos_communicator),
				      triangulation);

  {
    stokes_dof_handler.distribute_dofs (stokes_fe);
    DoFRenumbering::subdomain_wise (stokes_dof_handler);
    DoFRenumbering::component_wise (stokes_dof_handler, stokes_sub_blocks);

    stokes_constraints.clear ();
    DoFTools::make_hanging_node_constraints (stokes_dof_handler,
					     stokes_constraints);
    std::set<unsigned char> no_normal_flux_boundaries;
    no_normal_flux_boundaries.insert (0);
    VectorTools::compute_no_normal_flux_constraints (stokes_dof_handler, 0,
						     no_normal_flux_boundaries,
						     stokes_constraints);
    stokes_constraints.close ();
  }
  {
    temperature_dof_handler.distribute_dofs (temperature_fe);
    DoFRenumbering::subdomain_wise (temperature_dof_handler);

    temperature_constraints.clear ();
    DoFTools::make_hanging_node_constraints (temperature_dof_handler,
					     temperature_constraints);
    temperature_constraints.close ();
  }

  std::vector<unsigned int> stokes_dofs_per_block (2);
  DoFTools::count_dofs_per_block (stokes_dof_handler, stokes_dofs_per_block,
				  stokes_sub_blocks);

  const unsigned int n_u = stokes_dofs_per_block[0],
                     n_p = stokes_dofs_per_block[1],
		     n_T = temperature_dof_handler.n_dofs();

  pcout << "Number of active cells: "
	<< triangulation.n_active_cells()
	<< " (on "
	<< triangulation.n_levels()
	<< " levels)"
	<< std::endl
	<< "Number of degrees of freedom: "
	<< n_u + n_p + n_T
	<< " (" << n_u << '+' << n_p << '+'<< n_T <<')'
	<< std::endl
	<< std::endl;

  stokes_partitioner.clear();
  {
    std::vector<unsigned int> local_dofs (dim+1);
    DoFTools::
      count_dofs_with_subdomain_association (stokes_dof_handler,
					     Utilities::Trilinos::
					     get_this_mpi_process(trilinos_communicator),
					     local_dofs);
    const unsigned int
      n_local_velocities = std::accumulate (&local_dofs[0],
					    &local_dofs[dim],
					    0),
      n_local_pressures  = local_dofs[dim];

    stokes_partitioner.push_back (Epetra_Map(n_u, n_local_velocities,
					     0, trilinos_communicator));
    stokes_partitioner.push_back (Epetra_Map(n_p, n_local_pressures,
					     0, trilinos_communicator));
  }

  temperature_partitioner
    = Epetra_Map (n_T,
		  DoFTools::count_dofs_with_subdomain_association
		  (temperature_dof_handler,
		   Utilities::Trilinos::get_this_mpi_process(trilinos_communicator)),
		  0,
		  trilinos_communicator);

  if (Utilities::Trilinos::get_n_mpi_processes(trilinos_communicator) == 1)
    {
      Threads::TaskGroup<> tasks;
      tasks += Threads::new_task (&BoussinesqFlowProblem<dim>::setup_stokes_matrix,
				  *this);
      tasks += Threads::new_task (&BoussinesqFlowProblem<dim>::setup_stokes_preconditioner,
				  *this);
      tasks += Threads::new_task (&BoussinesqFlowProblem<dim>::setup_temperature_matrices,
				  *this);
      tasks.join_all ();
    }
  else
    {
      setup_stokes_matrix ();
      setup_stokes_preconditioner ();
      setup_temperature_matrices ();
    }

  stokes_solution.reinit (stokes_partitioner);
  old_stokes_solution.reinit (stokes_partitioner);
  stokes_rhs.reinit (stokes_partitioner);

  temperature_solution.reinit (temperature_partitioner);
  old_temperature_solution.reinit (temperature_partitioner);
  old_old_temperature_solution.reinit (temperature_partitioner);
  temperature_rhs.reinit (temperature_partitioner);

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
				     scratch.phi_p[i] * scratch.phi_p[j])
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
			  (Utilities::Trilinos::get_this_mpi_process(trilinos_communicator)),
			  stokes_dof_handler.begin_active()),
	 SubdomainFilter (IteratorFilters::SubdomainEqualTo
			  (Utilities::Trilinos::get_this_mpi_process(trilinos_communicator)),
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
  Amg_data.aggregation_threshold = 0.02;

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
				       - scratch.div_phi_u[i] * scratch.phi_p[j]
				       - scratch.phi_p[i] * scratch.div_phi_u[j])
				      * scratch.stokes_fe_values.JxW(q);

      const Point<dim> gravity = ( (dim == 2) ? (Point<dim> (0,1)) :
				   (Point<dim> (0,0,1)) );
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	data.local_rhs(i) += (EquationData::Rayleigh_number *
			      gravity * scratch.phi_u[i] * old_temperature)*
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
  pcout << "   Assembling..." << std::flush;

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
			  (Utilities::Trilinos::get_this_mpi_process(trilinos_communicator)),
			  stokes_dof_handler.begin_active()),
	 SubdomainFilter (IteratorFilters::SubdomainEqualTo
			  (Utilities::Trilinos::get_this_mpi_process(trilinos_communicator)),
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
  stokes_rhs.compress();

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
			  (Utilities::Trilinos::get_this_mpi_process(trilinos_communicator)),
			  temperature_dof_handler.begin_active()),
	 SubdomainFilter (IteratorFilters::SubdomainEqualTo
			  (Utilities::Trilinos::get_this_mpi_process(trilinos_communicator)),
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
				 // to the other assembly functions.
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

  EquationData::TemperatureRightHandSide<dim>  temperature_right_hand_side;

  const FEValuesExtractors::Vector velocities (0);

  data.local_rhs = 0;

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

  temperature_right_hand_side.value_list (scratch.temperature_fe_values.get_quadrature_points(),
					  scratch.gamma_values);

  scratch.stokes_fe_values[velocities].get_function_values (stokes_solution,
							    scratch.old_velocity_values);
  scratch.stokes_fe_values[velocities].get_function_values (old_stokes_solution,
							    scratch.old_old_velocity_values);

  const double nu
    = compute_viscosity (scratch.old_temperature_values,
			 scratch.old_old_temperature_values,
			 scratch.old_temperature_grads,
			 scratch.old_old_temperature_grads,
			 scratch.old_temperature_laplacians,
			 scratch.old_old_temperature_laplacians,
			 scratch.old_velocity_values,
			 scratch.old_old_velocity_values,
			 scratch.gamma_values,
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

      for (unsigned int i=0; i<dofs_per_cell; ++i)
	data.local_rhs(i) += (old_Ts * scratch.phi_T[i]
			      -
			      time_step *
			      extrapolated_u * ext_grad_T * scratch.phi_T[i]
			      -
			      time_step *
			      nu * ext_grad_T * scratch.grad_phi_T[i]
			      +
			      time_step *
			      scratch.gamma_values[q] * scratch.phi_T[i])
			     *
			     scratch.temperature_fe_values.JxW(q);
    }

  cell->get_dof_indices (data.local_dof_indices);
}


template <int dim>
void
BoussinesqFlowProblem<dim>::
copy_local_to_global_temperature_rhs (const Assembly::CopyData::TemperatureRHS<dim> &data)
{
  temperature_constraints.distribute_local_to_global (data.local_rhs,
						      data.local_dof_indices,
						      temperature_rhs);
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
			  (Utilities::Trilinos::get_this_mpi_process(trilinos_communicator)),
			  temperature_dof_handler.begin_active()),
	 SubdomainFilter (IteratorFilters::SubdomainEqualTo
			  (Utilities::Trilinos::get_this_mpi_process(trilinos_communicator)),
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

  temperature_rhs.compress();
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
  pcout << "   Solving..." << std::endl;

  {
    const LinearSolvers::BlockSchurPreconditioner<TrilinosWrappers::PreconditionAMG,
                                                  TrilinosWrappers::PreconditionILU>
      preconditioner (stokes_matrix, *Mp_preconditioner, *Amg_preconditioner);

    TrilinosWrappers::MPI::BlockVector
      distributed_stokes_solution (stokes_partitioner);
    distributed_stokes_solution = stokes_solution;

    const unsigned int
      start = (distributed_stokes_solution.block(0).size() +
	       distributed_stokes_solution.block(1).local_range().first),
      end   = (distributed_stokes_solution.block(0).size() +
	       distributed_stokes_solution.block(1).local_range().second);
    for (unsigned int i=start; i<end; ++i)
      if (stokes_constraints.is_constrained (i))
	distributed_stokes_solution(i) = 0;


    SolverControl solver_control (stokes_matrix.m(), 1e-6*stokes_rhs.l2_norm());
    SolverBicgstab<TrilinosWrappers::MPI::BlockVector> 
      bicgstab (solver_control, false);

    bicgstab.solve(stokes_matrix, distributed_stokes_solution, stokes_rhs,
		   preconditioner);

    stokes_solution = distributed_stokes_solution;

    pcout << "   "
	  << solver_control.last_step()
	  << " BiCGStab iterations for Stokes subsystem."
	  << std::endl;

    stokes_constraints.distribute (stokes_solution);
  }
  computing_timer.exit_section();


  computing_timer.enter_section ("   Assemble temperature rhs");

  old_time_step = time_step;
  const double maximal_velocity = get_maximal_velocity();
  time_step = 1./(1.6*dim*std::sqrt(1.*dim)) /
	      temperature_degree *
	      GridTools::minimal_cell_diameter(triangulation) /
              std::max (maximal_velocity, 0.01);

  pcout << "   " << "Time step: " << time_step
	<< std::endl;

  temperature_solution = old_temperature_solution;

  assemble_temperature_system (maximal_velocity);

  computing_timer.exit_section ();

  computing_timer.enter_section ("   Solve temperature system");

  {
    SolverControl solver_control (temperature_matrix.m(),
				  1e-8*temperature_rhs.l2_norm());
    SolverCG<TrilinosWrappers::MPI::Vector>   cg (solver_control);

    TrilinosWrappers::MPI::Vector
      distributed_temperature_solution (temperature_partitioner);
    distributed_temperature_solution = temperature_solution;

    cg.solve (temperature_matrix, distributed_temperature_solution,
	      temperature_rhs, *T_preconditioner);

    temperature_solution = distributed_temperature_solution;
    temperature_constraints.distribute (temperature_solution);

    pcout << "   "
	  << solver_control.last_step()
	  << " CG iterations for temperature" << std::endl;
    computing_timer.exit_section();

    double min_temperature = temperature_solution(0),
	   max_temperature = temperature_solution(0);
    for (unsigned int i=1; i<temperature_solution.size(); ++i)
      {
	min_temperature = std::min<double> (min_temperature,
					    temperature_solution(i));
	max_temperature = std::max<double> (max_temperature,
					    temperature_solution(i));
      }

    pcout << "   Temperature range: "
	  << min_temperature << ' ' << max_temperature
	  << std::endl;
  }
}



				 // @sect4{BoussinesqFlowProblem::output_results}

				 // This function has remained completely
				 // unchanged compared to step-31 (with the
				 // exception that we make sure that only a
				 // single processor actually does some work
				 // here), so everything should be clear here:
template <int dim>
void BoussinesqFlowProblem<dim>::output_results ()
{
  if (timestep_number % 10 != 0)
    return;

  computing_timer.enter_section ("Postprocessing");
  
  if (Utilities::Trilinos::get_this_mpi_process(trilinos_communicator) == 0)
    {

      const FESystem<dim> joint_fe (stokes_fe, 1,
				    temperature_fe, 1);
      DoFHandler<dim> joint_dof_handler (triangulation);
      joint_dof_handler.distribute_dofs (joint_fe);
      Assert (joint_dof_handler.n_dofs() ==
	      stokes_dof_handler.n_dofs() + temperature_dof_handler.n_dofs(),
	      ExcInternalError());

      Vector<double> joint_solution (joint_dof_handler.n_dofs());
      
      {
	std::vector<unsigned int> local_joint_dof_indices (joint_fe.dofs_per_cell);
	std::vector<unsigned int> local_stokes_dof_indices (stokes_fe.dofs_per_cell);
	std::vector<unsigned int> local_temperature_dof_indices (temperature_fe.dofs_per_cell);

	typename DoFHandler<dim>::active_cell_iterator
	  joint_cell       = joint_dof_handler.begin_active(),
	  joint_endc       = joint_dof_handler.end(),
	  stokes_cell      = stokes_dof_handler.begin_active(),
	  temperature_cell = temperature_dof_handler.begin_active();
	for (; joint_cell!=joint_endc; ++joint_cell, ++stokes_cell, ++temperature_cell)
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
		    = stokes_solution(local_stokes_dof_indices[joint_fe.system_to_base_index(i).second]);
		}
	      else
		{
		  Assert (joint_fe.system_to_base_index(i).first.first == 1,
			  ExcInternalError());
		  Assert (joint_fe.system_to_base_index(i).second
			  <
			  local_stokes_dof_indices.size(),
			  ExcInternalError());
		  joint_solution(local_joint_dof_indices[i])
		    = temperature_solution(local_temperature_dof_indices[joint_fe.system_to_base_index(i).second]);
		}
	  }
      }

      std::vector<std::string> joint_solution_names (dim, "velocity");
      joint_solution_names.push_back ("p");
      joint_solution_names.push_back ("T");

      DataOut<dim> data_out;

      data_out.attach_dof_handler (joint_dof_handler);

      std::vector<DataComponentInterpretation::DataComponentInterpretation>
	data_component_interpretation
	(dim+2, DataComponentInterpretation::component_is_scalar);
      for (unsigned int i=0; i<dim; ++i)
	data_component_interpretation[i]
	  = DataComponentInterpretation::component_is_part_of_vector;

      data_out.add_data_vector (joint_solution, joint_solution_names,
				DataOut<dim>::type_dof_data,
				data_component_interpretation);
      data_out.build_patches (std::min(stokes_degree, temperature_degree));

      std::ostringstream filename;
      filename << "solution-" << Utilities::int_to_string(timestep_number, 4) << ".vtk";

      std::ofstream output (filename.str().c_str());
      data_out.write_vtk (output);
    }

  computing_timer.exit_section ();
}



				 // @sect4{BoussinesqFlowProblem::refine_mesh}

				 // Nothing new here, either. Since the
				 // <code>setup_dofs</code> function that we
				 // call in the middle has its own timer
				 // section, we split timing this function
				 // into two sections. It will also allow us
				 // to easily identify which of the two is
				 // more expensive.
template <int dim>
void BoussinesqFlowProblem<dim>::refine_mesh (const unsigned int max_grid_level)
{
  computing_timer.enter_section ("Refine mesh structure, part 1");
  Vector<float> estimated_error_per_cell (triangulation.n_active_cells());

  KellyErrorEstimator<dim>::estimate (temperature_dof_handler,
				      QGauss<dim-1>(temperature_degree+1),
				      typename FunctionMap<dim>::type(),
				      temperature_solution,
				      estimated_error_per_cell);

  GridRefinement::refine_and_coarsen_fixed_fraction (triangulation,
						     estimated_error_per_cell,
						     0.8, 0.1);
  if (triangulation.n_levels() > max_grid_level)
    for (typename Triangulation<dim>::active_cell_iterator
	   cell = triangulation.begin_active(max_grid_level);
	 cell != triangulation.end(); ++cell)
      cell->clear_refine_flag ();

  std::vector<TrilinosWrappers::Vector> x_temperature (2);
  x_temperature[0] = temperature_solution;
  x_temperature[1] = old_temperature_solution;
  TrilinosWrappers::BlockVector x_stokes = stokes_solution;

  SolutionTransfer<dim,TrilinosWrappers::Vector>
    temperature_trans(temperature_dof_handler);
  SolutionTransfer<dim,TrilinosWrappers::BlockVector>
    stokes_trans(stokes_dof_handler);

  triangulation.prepare_coarsening_and_refinement();
  temperature_trans.prepare_for_coarsening_and_refinement(x_temperature);
  stokes_trans.prepare_for_coarsening_and_refinement(x_stokes);

  triangulation.execute_coarsening_and_refinement ();
  computing_timer.exit_section();

  setup_dofs ();

  computing_timer.enter_section ("Refine mesh structure, part 2");

  std::vector<TrilinosWrappers::Vector> tmp (2);
  tmp[0].reinit (temperature_solution);
  tmp[1].reinit (temperature_solution);
  temperature_trans.interpolate(x_temperature, tmp);

  temperature_solution = tmp[0];
  old_temperature_solution = tmp[1];

  TrilinosWrappers::BlockVector x_stokes_new = stokes_solution;
  stokes_trans.interpolate (x_stokes, x_stokes_new);
  stokes_solution = x_stokes_new;

  rebuild_stokes_matrix              = true;
  rebuild_stokes_preconditioner      = true;
  rebuild_temperature_matrices       = true;
  rebuild_temperature_preconditioner = true;

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
void BoussinesqFlowProblem<dim>::run ()
{
  const unsigned int initial_refinement = (dim == 2 ? 4 : 2);
  const unsigned int n_pre_refinement_steps = (dim == 2 ? 4 : 3);

				 /* Data for the shell problem. */
  /*
  GridGenerator::half_hyper_shell (triangulation,
				   Point<dim>(), 0.5, 1.0);

  static HyperShellBoundary<dim> boundary;
  triangulation.set_boundary (0, boundary);
  */

				 /* Data for the cube problem. */
  GridGenerator::hyper_cube (triangulation);
  global_Omega_diameter = GridTools::diameter (triangulation);

  triangulation.refine_global (initial_refinement);

  setup_dofs();

  unsigned int       pre_refinement_step    = 0;

  start_time_iteration:

  project_temperature_field ();

  timestep_number           = 0;
  time_step = old_time_step = 0;

  double time = 0;

  do
    {
      pcout << "Timestep " << timestep_number
	    << ":  t=" << time
	    << std::endl;

      assemble_stokes_system ();
      build_stokes_preconditioner ();
      assemble_temperature_matrix ();

      solve ();

      output_results ();

      pcout << std::endl;

      if ((timestep_number == 0) &&
	  (pre_refinement_step < n_pre_refinement_steps))
	{
	  refine_mesh (initial_refinement + n_pre_refinement_steps);
	  ++pre_refinement_step;
	  goto start_time_iteration;
	}
      else
	if ((timestep_number > 0) && (timestep_number % 5 == 0))
	  refine_mesh (initial_refinement + n_pre_refinement_steps);

      time += time_step;
      ++timestep_number;

      old_stokes_solution          = stokes_solution;
      old_old_temperature_solution = old_temperature_solution;
      old_temperature_solution     = temperature_solution;
    }
  while (time <= 100);
}



				 // @sect3{The <code>main</code> function}

				 // This is copied verbatim from step-31:
int main (int argc, char *argv[])
{
  try
    {
      deallog.depth_console (0);

      Utilities::System::MPI_InitFinalize mpi_initialization(argc, argv);

      BoussinesqFlowProblem<2> flow_problem;
      flow_problem.run ();
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
