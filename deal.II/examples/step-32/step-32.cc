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

#include <Epetra_Map.h>

				 // Time measurements. 
#include <base/timer.h>

#include <fstream>
#include <iostream>
#include <sstream>


				 // Next, we import all deal.II
				 // names into global namespace
using namespace dealii;

				 // @sect3{Equation data}


				 // @sect4{Boundary values}
namespace EquationData
{
				   // define viscosity
  const double eta = 1;
  const double kappa = 1e-6;
  const double Rayleigh_number = 10;


				   // @sect4{Initial values}
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
    //return (p.norm() < 0.55+0.02*std::sin(p[0]*20) ? 1 : 0);
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



				   // @sect4{Right hand side}
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
    //    return 0;
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
    void output_results () const;
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
    std_cxx1x::shared_ptr<TrilinosWrappers::PreconditionILU>  Mp_preconditioner;
    std_cxx1x::shared_ptr<TrilinosWrappers::PreconditionILU>  T_preconditioner;

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
    
    void copy_local_to_global_stokes_preconditioner (const Assembly::CopyData::StokesPreconditioner<dim> &data);


    void
    local_assemble_stokes_system (const typename DoFHandler<dim>::active_cell_iterator &cell,
				  Assembly::Scratch::StokesSystem<dim>  &scratch,
				  Assembly::CopyData::StokesSystem<dim> &data);
    
    void copy_local_to_global_stokes_system (const Assembly::CopyData::StokesSystem<dim> &data);


    void
    local_assemble_temperature_matrix (const typename DoFHandler<dim>::active_cell_iterator &cell,
				       Assembly::Scratch::TemperatureMatrix<dim>  &scratch,
				       Assembly::CopyData::TemperatureMatrix<dim> &data);
    
    void copy_local_to_global_temperature_matrix (const Assembly::CopyData::TemperatureMatrix<dim> &data);



    void
    local_assemble_temperature_rhs (const std::pair<double,double> global_T_range,
				    const double                   global_max_velocity,
				    const typename DoFHandler<dim>::active_cell_iterator &cell,
				    Assembly::Scratch::TemperatureRHS<dim> &scratch,
				    Assembly::CopyData::TemperatureRHS<dim> &data);
    
    void copy_local_to_global_temperature_rhs (const Assembly::CopyData::TemperatureRHS<dim> &data);
};


				 // @sect3{BoussinesqFlowProblem class implementation}

				 // @sect4{BoussinesqFlowProblem::BoussinesqFlowProblem}
template <int dim>
BoussinesqFlowProblem<dim>::BoussinesqFlowProblem ()
                :
                trilinos_communicator (Utilities::Trilinos::comm_world()),
		pcout (std::cout, Utilities::Trilinos::get_this_mpi_process(trilinos_communicator)==0),

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



				 // @sect4{BoussinesqFlowProblem::get_maximal_velocity}
template <int dim>
double BoussinesqFlowProblem<dim>::get_maximal_velocity () const
{
  const QIterated<dim> quadrature_formula (QTrapez<1>(),
					   stokes_degree+1);
  const unsigned int n_q_points = quadrature_formula.size();

  FEValues<dim> fe_values (stokes_fe, quadrature_formula, update_values);
  std::vector<Tensor<1,dim> > velocity_values(n_q_points);

  const FEValuesExtractors::Vector velocities (0);
  
  double max_local_velocity = 0, max_velocity = 0;

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

  trilinos_communicator.MaxAll(&max_local_velocity, &max_velocity, 1);

  return max_velocity;
}




				 // @sect4{BoussinesqFlowProblem::get_extrapolated_temperature_range}
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
				   /*std::cout << "Processor " << trilinos_communicator.MyPID() 
				     << " stokes(0,0) rows: " 
				     << stokes_matrix.block(0,0).matrix->NumMyRows()
				     << ", nnz: " 
				     << stokes_matrix.block(0,0).matrix->NumMyNonzeros() 
				     << std::endl;*/

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




				 // @sect4{BoussinesqFlowProblem::setup_dofs}
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
    unsigned int n_local_velocities = 0;
    for (unsigned int c=0; c<dim; ++c)
      n_local_velocities += local_dofs[c];

    const unsigned int n_local_pressures = local_dofs[dim];

    Epetra_Map map_u(n_u, n_local_velocities, 0, trilinos_communicator);
    stokes_partitioner.push_back (map_u);
    Epetra_Map map_p(n_p, n_local_pressures, 0, trilinos_communicator);
    stokes_partitioner.push_back (map_p);
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
				     scalar_product (scratch.grad_phi_u[i], scratch.grad_phi_u[j])
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
				       scratch.grads_phi_u[i] * scratch.grads_phi_u[j]
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



// @sect4{BoussinesqFlowProblem::assemble_stokes_system}
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






				 // @sect4{BoussinesqFlowProblem::assemble_temperature_system}
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
      T_preconditioner =  std_cxx1x::shared_ptr<TrilinosWrappers::PreconditionILU>
                                   (new TrilinosWrappers::PreconditionILU());
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



				 // @sect4{BoussinesqFlowProblem::project_temperature_field}
				 // Manually project the initial
				 // conditions for the temperature in
				 // parallel instead of doing that
				 // completely on each processor. The
				 // temperature mass matrix is already
				 // available, and we need just to
				 // compute a right hand side in
				 // parallel, do a cg solve and
				 // distribute the hanging node
				 // constraints.
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

  TrilinosWrappers::MPI::Vector rhs (temperature_mass_matrix.row_partitioner()), 
    sol (temperature_mass_matrix.row_partitioner());
      
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

  TrilinosWrappers::PreconditionILU preconditioner_mass;
  preconditioner_mass.initialize(temperature_mass_matrix);

  cg.solve (temperature_mass_matrix, sol, rhs, preconditioner_mass);
  
  old_temperature_solution = sol;
  temperature_constraints.distribute (old_temperature_solution);
}

				 // @sect4{BoussinesqFlowProblem::solve}
template <int dim>
void BoussinesqFlowProblem<dim>::solve ()
{
  computing_timer.enter_section ("   Solve Stokes system");
  pcout << "   Solving..." << std::endl;

  {
    const LinearSolvers::BlockSchurPreconditioner<TrilinosWrappers::PreconditionAMG,
                                                  TrilinosWrappers::PreconditionILU>
      preconditioner (stokes_matrix, *Mp_preconditioner, *Amg_preconditioner);

    SolverControl solver_control (stokes_matrix.m(),
				  1e-6*stokes_rhs.l2_norm());

    SolverBicgstab<TrilinosWrappers::MPI::BlockVector>
      bicgstab (solver_control, false);

    TrilinosWrappers::MPI::BlockVector
      distributed_stokes_solution (stokes_partitioner);
    distributed_stokes_solution = stokes_solution;

				   // now treat the hanging nodes correctly.
    const unsigned int start = 
      distributed_stokes_solution.block(1).local_range().first + 
      distributed_stokes_solution.block(0).size();
    const unsigned int end = 
      distributed_stokes_solution.block(1).local_range().second + 
      distributed_stokes_solution.block(0).size();
    for (unsigned int i=start; i<end; ++i)
      if (stokes_constraints.is_constrained (i))
	distributed_stokes_solution(i) = 0;
    
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
template <int dim>
void BoussinesqFlowProblem<dim>::output_results ()  const
{
  if (timestep_number % 10 != 0)
    return;

  const FESystem<dim> joint_fe (stokes_fe, 1,
				temperature_fe, 1);
  DoFHandler<dim> joint_dof_handler (triangulation);
  joint_dof_handler.distribute_dofs (joint_fe);
  Assert (joint_dof_handler.n_dofs() ==
	  stokes_dof_handler.n_dofs() + temperature_dof_handler.n_dofs(),
	  ExcInternalError());

  Vector<double> joint_solution (joint_dof_handler.n_dofs());

  if (Utilities::Trilinos::get_this_mpi_process(trilinos_communicator) == 0)
    {

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
}



				 // @sect4{BoussinesqFlowProblem::refine_mesh}
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
template <int dim>
void BoussinesqFlowProblem<dim>::run ()
{
  const unsigned int initial_refinement = (dim == 2 ? 4 : 2);
  const unsigned int n_pre_refinement_steps = (dim == 2 ? 4 : 3);

  //GridGenerator::half_hyper_shell (triangulation,
  //				   Point<dim>(), 0.5, 1.0);

  //static HyperShellBoundary<dim> boundary;
  //triangulation.set_boundary (0, boundary);
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
