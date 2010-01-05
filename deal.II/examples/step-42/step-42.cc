/* $Id$ */
/* Author: Wolfgang Bangerth, Texas A&M University, 2008 */

/*    $Id$       */
/*                                                                */
/*    Copyright (C) 2008, 2009, 2010 by the deal.II authors */
/*                                                                */
/*    This file is subject to QPL and may not be  distributed     */
/*    without copyright and license information. Please refer     */
/*    to the file deal.II/doc/license.html for the  text  and     */
/*    further information on this license.                        */


				 // @sect3{Include files}
                        
				 // As usual, we start by including 
				 // some well-known files:
#include <base/quadrature_lib.h>
#include <base/logstream.h>
#include <base/function.h>
#include <base/utilities.h>

#include <lac/block_vector.h>
#include <lac/full_matrix.h>
#include <lac/block_sparse_matrix.h>
#include <lac/solver_cg.h>
#include <lac/solver_gmres.h>
#include <lac/precondition.h>
#include <lac/constraint_matrix.h>

#include <grid/tria.h>
#include <grid/grid_generator.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
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

				
#include <lac/sparse_direct.h>

#include <lac/sparse_ilu.h>

#include <multigrid/multigrid.h>
#include <multigrid/mg_dof_handler.h>
#include <multigrid/mg_dof_accessor.h>
#include <multigrid/mg_transfer.h>
#include <multigrid/mg_tools.h>
#include <multigrid/mg_coarse.h>
#include <multigrid/mg_smoother.h>
#include <multigrid/mg_matrix.h>
				
#include <fstream>
#include <sstream>

				
using namespace dealii;

			
template <int dim>
struct InnerPreconditioner;

			
template <>
struct InnerPreconditioner<2> 
{
    typedef SparseDirectUMFPACK type;
};

				
template <>
struct InnerPreconditioner<3> 
{
    typedef SparseILU<double> type;
};


template <int dim>
class StokesProblem 
{
  public:
    StokesProblem (const unsigned int degree);
    void run ();
    
  private:
    void setup_dofs ();
    void assemble_system ();
    void assemble_multigrid ();
    void solve_schur ();
    void solve_block ();
    void solve_mg ();

    void find_dofs_on_lower_level (std::vector<std::vector<bool> > &lower_dofs, 
        std::vector<std::vector<bool> > &boundary_dofs);
    void output_results (const unsigned int refinement_cycle) const;
    void refine_mesh ();
    
    const unsigned int   degree;
    
    Triangulation<dim>   triangulation;
    FESystem<dim>        fe;
    MGDoFHandler<dim>    dof_handler;

    ConstraintMatrix     constraints;
    
    BlockSparsityPattern      sparsity_pattern;
    BlockSparseMatrix<double> system_matrix;

    BlockVector<double> solution;
    BlockVector<double> system_rhs;

    MGLevelObject<ConstraintMatrix>           mg_constraints;
    MGLevelObject<BlockSparsityPattern>       mg_sparsity;
    MGLevelObject<BlockSparseMatrix<double> > mg_matrices;

    MGLevelObject<BlockSparseMatrix<double> > mg_interface_matrices;
    std::vector<std::vector<unsigned int> >   mg_dofs_per_component;
				  
    std_cxx1x::shared_ptr<typename InnerPreconditioner<dim>::type> A_preconditioner;
    std::vector<std_cxx1x::shared_ptr<typename InnerPreconditioner<dim>::type> > mg_A_preconditioner;
};

			
template <int dim>
class BoundaryValues : public Function<dim> 
{
  public:
    BoundaryValues () : Function<dim>(dim+1) {}
    
    virtual double value (const Point<dim>   &p,
                          const unsigned int  component = 0) const;

    virtual void vector_value (const Point<dim> &p, 
                               Vector<double>   &value) const;
};


template <int dim>
double
BoundaryValues<dim>::value (const Point<dim>  &p,
			    const unsigned int component) const 
{
  Assert (component < this->n_components,
	  ExcIndexRange (component, 0, this->n_components));
  
   if (component == 0 && p[0] == 0)
    return (dim == 2 ? - p[1]*(p[1]-1.) : p[1]*(p[1]-1.) * p[2]*(p[2]-1.));
  return 0;
}


template <int dim>
void
BoundaryValues<dim>::vector_value (const Point<dim> &p,
				   Vector<double>   &values) const 
{
  for (unsigned int c=0; c<this->n_components; ++c)
    values(c) = BoundaryValues<dim>::value (p, c);
}



			
template <int dim>
class RightHandSide : public Function<dim> 
{
  public:
    RightHandSide () : Function<dim>(dim+1) {}
    
    virtual double value (const Point<dim>   &p,
                          const unsigned int  component = 0) const;

    virtual void vector_value (const Point<dim> &p, 
                               Vector<double>   &value) const;
    
};


template <int dim>
double
RightHandSide<dim>::value (const Point<dim>  &/*p*/,
                           const unsigned int /*component*/) const 
{
  return 0;
}


template <int dim>
void
RightHandSide<dim>::vector_value (const Point<dim> &p,
                                  Vector<double>   &values) const 
{
  for (unsigned int c=0; c<this->n_components; ++c)
    values(c) = RightHandSide<dim>::value (p, c);
}


			
				
template <class Matrix, class Preconditioner>
class InverseMatrix : public Subscriptor
{
  public:
    InverseMatrix (const Matrix         &m,
                   const Preconditioner &preconditioner);

    void vmult (Vector<double>       &dst,
                const Vector<double> &src) const;

  private:
    const SmartPointer<const Matrix> matrix;
    const SmartPointer<const Preconditioner> preconditioner;
};


template <class Matrix, class Preconditioner>
InverseMatrix<Matrix,Preconditioner>::InverseMatrix (const Matrix &m,
						     const Preconditioner &preconditioner)
		:
		matrix (&m),
		preconditioner (&preconditioner)
{}


template <class Matrix, class Preconditioner>
void InverseMatrix<Matrix,Preconditioner>::vmult (Vector<double>       &dst,
						  const Vector<double> &src) const
{
  SolverControl solver_control (src.size(), 1.0e-12);
  SolverCG<>    cg (solver_control);

  dst = 0;

  cg.solve (*matrix, dst, src, *preconditioner);

  std::cout << "  Preconditioner CG steps: " << solver_control.last_step() << std::endl;
}

template <class PreconditionerA, class PreconditionerMp>
class BlockSchurPreconditioner : public Subscriptor
{
  public:
    BlockSchurPreconditioner (const BlockSparseMatrix<double>         &S,
          const InverseMatrix<SparseMatrix<double>,PreconditionerMp>  &Mpinv,
          const PreconditionerA &Apreconditioner);

  void vmult (BlockVector<double>       &dst,
              const BlockVector<double> &src) const;

  private:
    const SmartPointer<const BlockSparseMatrix<double> > system_matrix;
    const SmartPointer<const InverseMatrix<SparseMatrix<double>, 
                       PreconditionerMp > > m_inverse;
    const PreconditionerA &a_preconditioner;
    
    mutable Vector<double> tmp;

};

template <class PreconditionerA, class PreconditionerMp>
BlockSchurPreconditioner<PreconditionerA, PreconditionerMp>::BlockSchurPreconditioner(
          const BlockSparseMatrix<double>                            &S,
          const InverseMatrix<SparseMatrix<double>,PreconditionerMp> &Mpinv,
          const PreconditionerA &Apreconditioner
          )
                :
                system_matrix           (&S),
                m_inverse               (&Mpinv),
                a_preconditioner        (Apreconditioner),
                tmp                     (S.block(1,1).m())
{}

        // Now the interesting function, the multiplication of
        // the preconditioner with a BlockVector. 
template <class PreconditionerA, class PreconditionerMp>
void BlockSchurPreconditioner<PreconditionerA, PreconditionerMp>::vmult (
                                     BlockVector<double>       &dst,
                                     const BlockVector<double> &src) const
{
        // Form u_new = A^{-1} u
  a_preconditioner.vmult (dst.block(0), src.block(0));
        // Form tmp = - B u_new + p 
        // (<code>SparseMatrix::residual</code>
        // does precisely this)
  system_matrix->block(1,0).residual(tmp, dst.block(0), src.block(1));
        // Change sign in tmp
  tmp *= -1;
        // Multiply by approximate Schur complement 
        // (i.e. a pressure mass matrix)
  m_inverse->vmult (dst.block(1), tmp);
}


template <class Preconditioner>
class SchurComplement : public Subscriptor
{
  public:
    SchurComplement (const BlockSparseMatrix<double> &system_matrix,
		     const InverseMatrix<SparseMatrix<double>, Preconditioner> &A_inverse);

    void vmult (Vector<double>       &dst,
		const Vector<double> &src) const;

  private:
    const SmartPointer<const BlockSparseMatrix<double> > system_matrix;
    const SmartPointer<const InverseMatrix<SparseMatrix<double>, Preconditioner> > A_inverse;
    
    mutable Vector<double> tmp1, tmp2;
};



template <class Preconditioner>
SchurComplement<Preconditioner>::
SchurComplement (const BlockSparseMatrix<double> &system_matrix,
		 const InverseMatrix<SparseMatrix<double>,Preconditioner> &A_inverse)
		:
		system_matrix (&system_matrix),
		A_inverse (&A_inverse),
		tmp1 (system_matrix.block(0,0).m()),
		tmp2 (system_matrix.block(0,0).m())
{}


template <class Preconditioner>
void SchurComplement<Preconditioner>::vmult (Vector<double>       &dst,
					     const Vector<double> &src) const
{
  system_matrix->block(0,1).vmult (tmp1, src);
  A_inverse->vmult (tmp2, tmp1);
  system_matrix->block(1,0).vmult (dst, tmp2);
}


			
template <int dim>
StokesProblem<dim>::StokesProblem (const unsigned int degree)
                :
                degree (degree),
                triangulation (Triangulation<dim>::maximum_smoothing),
                fe (FE_Q<dim>(degree+1), dim,
                    FE_Q<dim>(degree), 1),
                dof_handler (triangulation)
{}


				
			
				
template <int dim>
void StokesProblem<dim>::setup_dofs ()
{
  A_preconditioner.reset ();
  mg_A_preconditioner.resize (0);
  system_matrix.clear ();
  
  dof_handler.distribute_dofs (fe);  
//  DoFRenumbering::Cuthill_McKee (dof_handler);

  std::vector<unsigned int> block_component (dim+1,0);
  block_component[dim] = 1;
  DoFRenumbering::component_wise (dof_handler, block_component);

				
  {
    constraints.clear ();
    std::vector<bool> component_mask (dim+1, true);
    component_mask[dim] = false;
    VectorTools::interpolate_boundary_values (dof_handler,
					      0,
					      BoundaryValues<dim>(),
					      constraints,
					      component_mask);
    DoFTools::make_hanging_node_constraints (dof_handler,
					     constraints);
  }

  constraints.close ();

				
  std::vector<unsigned int> dofs_per_block (2);
  DoFTools::count_dofs_per_block (dof_handler, dofs_per_block, block_component);  
  const unsigned int n_u = dofs_per_block[0],
                     n_p = dofs_per_block[1];

  std::cout << "   Number of active cells: "
            << triangulation.n_active_cells()
            << std::endl
            << "   Number of degrees of freedom: "
            << dof_handler.n_dofs()
            << " (" << n_u << '+' << n_p << ')'
            << std::endl;
      
				
				 
				 
  {
    BlockCompressedSimpleSparsityPattern csp (2,2);

    csp.block(0,0).reinit (n_u, n_u);
    csp.block(1,0).reinit (n_p, n_u);
    csp.block(0,1).reinit (n_u, n_p);
    csp.block(1,1).reinit (n_p, n_p);
  
    csp.collect_sizes();    
  
    DoFTools::make_sparsity_pattern (
        static_cast<const DoFHandler<dim>&>(dof_handler), 
        csp, constraints, false);
    sparsity_pattern.copy_from (csp);
  }
  
				 
  system_matrix.reinit (sparsity_pattern);
                                   
  solution.reinit (2);
  solution.block(0).reinit (n_u);
  solution.block(1).reinit (n_p);
  solution.collect_sizes ();
  
  system_rhs.reinit (2);
  system_rhs.block(0).reinit (n_u);
  system_rhs.block(1).reinit (n_p);
  system_rhs.collect_sizes ();

  //now setup stuff for mg
  const unsigned int nlevels = triangulation.n_levels();

  mg_matrices.resize(0, nlevels-1);
  mg_matrices.clear ();
  mg_interface_matrices.resize(0, nlevels-1);
  mg_interface_matrices.clear ();
  mg_sparsity.resize(0, nlevels-1);

  MGTools::count_dofs_per_block (dof_handler, mg_dofs_per_component);

  for (unsigned int level=0;level<nlevels;++level)
  {
    BlockCompressedSparsityPattern bscp (mg_dofs_per_component[level], 
        mg_dofs_per_component[level]);
    MGTools::make_sparsity_pattern(dof_handler, bscp, level);
    mg_sparsity[level].copy_from (bscp);
    mg_matrices[level].reinit (mg_sparsity[level]);
    if(level>0)
      mg_interface_matrices[level].reinit (mg_sparsity[level]);
  }
}


			
template <int dim>
void StokesProblem<dim>::assemble_system () 
{
  system_matrix=0;
  system_rhs=0;
  
  QGauss<dim>   quadrature_formula(degree+2);

  FEValues<dim> fe_values (fe, quadrature_formula,
                           update_values    |
                           update_quadrature_points  |
                           update_JxW_values |
                           update_gradients);
  
  const unsigned int   dofs_per_cell   = fe.dofs_per_cell;
  
  const unsigned int   n_q_points      = quadrature_formula.size();

  FullMatrix<double>   local_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>       local_rhs (dofs_per_cell);

  std::vector<unsigned int> local_dof_indices (dofs_per_cell);
  
  const RightHandSide<dim>          right_hand_side;
  std::vector<Vector<double> >      rhs_values (n_q_points,
                                                Vector<double>(dim+1));

				 
  const FEValuesExtractors::Vector velocities (0);
  const FEValuesExtractors::Scalar pressure (dim);

				 
				  
  std::vector<Tensor<2,dim> >          phi_grads_u (dofs_per_cell);
  std::vector<double>                  div_phi_u   (dofs_per_cell);
  std::vector<double>                  phi_p       (dofs_per_cell);
				   
  typename MGDoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    { 
      fe_values.reinit (cell);
      local_matrix = 0;
      local_rhs = 0;
      
      right_hand_side.vector_value_list(fe_values.get_quadrature_points(),
                                        rhs_values);
      
      for (unsigned int q=0; q<n_q_points; ++q)
	{
	  for (unsigned int k=0; k<dofs_per_cell; ++k)
	    {
	      phi_grads_u[k] = fe_values[velocities].gradient (k, q);
	      div_phi_u[k]   = fe_values[velocities].divergence (k, q);
	      phi_p[k]       = fe_values[pressure].value (k, q);
	    }

	  for (unsigned int i=0; i<dofs_per_cell; ++i)
	    {
	      for (unsigned int j=0; j<dofs_per_cell; ++j)
		{
		  local_matrix(i,j) += (scalar_product(phi_grads_u[i], phi_grads_u[j])
					- div_phi_u[i] * phi_p[j]
					- phi_p[i] * div_phi_u[j]
					+ phi_p[i] * phi_p[j])
				       * fe_values.JxW(q);     

		}

	      const unsigned int component_i =
		fe.system_to_component_index(i).first;
	      local_rhs(i) += fe_values.shape_value(i,q) * 
			      rhs_values[q](component_i) *
			      fe_values.JxW(q);
	    }
	}

				      


      cell->get_dof_indices (local_dof_indices);
      constraints.distribute_local_to_global (local_matrix, local_rhs,
					      local_dof_indices, 
					      system_matrix, system_rhs);
    }
  
				
  std::cout << "   Computing preconditioner..." << std::endl << std::flush;
      
  A_preconditioner
    = std_cxx1x::shared_ptr<typename InnerPreconditioner<dim>::type>(new typename InnerPreconditioner<dim>::type());
  A_preconditioner->initialize (system_matrix.block(0,0),
				typename InnerPreconditioner<dim>::type::AdditionalData());

}


template <int dim>
void StokesProblem<dim>::assemble_multigrid () 
{
  QGauss<dim>   quadrature_formula(degree+2);
  FEValues<dim> fe_values (fe, quadrature_formula,
                           update_values    |
                           update_quadrature_points  |
                           update_JxW_values |
                           update_gradients);
  
  const unsigned int   dofs_per_cell   = fe.dofs_per_cell;
  const unsigned int   n_q_points      = quadrature_formula.size();

  FullMatrix<double>   local_matrix (dofs_per_cell, dofs_per_cell);

  std::vector<unsigned int> local_dof_indices (dofs_per_cell);
				 
  const FEValuesExtractors::Vector velocities (0);
  const FEValuesExtractors::Scalar pressure (dim);
				 
				  
  std::vector<Tensor<2,dim> >          phi_grads_u (dofs_per_cell);
  std::vector<double>                  div_phi_u   (dofs_per_cell);
  std::vector<double>                  phi_p       (dofs_per_cell);
				   
  std::vector<std::vector<bool> > interface_dofs;
  std::vector<std::vector<bool> > boundary_interface_dofs;
  for (unsigned int level = 0; level<triangulation.n_levels(); ++level)
    {
      std::vector<bool> tmp (dof_handler.n_dofs(level));
      interface_dofs.push_back (tmp);
      boundary_interface_dofs.push_back (tmp);
    }
  MGTools::extract_inner_interface_dofs (dof_handler,
					 interface_dofs, boundary_interface_dofs);

  typename FunctionMap<dim>::type      dirichlet_boundary;
  BoundaryValues<dim>                  dirichlet_bc;
  dirichlet_boundary[0] =             &dirichlet_bc;

  std::vector<std::set<unsigned int> > boundary_indices(triangulation.n_levels());
  std::vector<bool> component_mask (dim+1, true);
  component_mask[dim] = false;
  MGTools::make_boundary_list (dof_handler, dirichlet_boundary,
			       boundary_indices, component_mask);

  std::vector<ConstraintMatrix> boundary_constraints (triangulation.n_levels());
  std::vector<ConstraintMatrix> boundary_interface_constraints (triangulation.n_levels());
  for (unsigned int level=0; level<triangulation.n_levels(); ++level)
    {
      boundary_constraints[level].add_lines (interface_dofs[level]);
      boundary_constraints[level].add_lines (boundary_indices[level]);
      boundary_constraints[level].close ();

      boundary_interface_constraints[level]
	.add_lines (boundary_interface_dofs[level]);
      boundary_interface_constraints[level].close ();
    }

  typename MGDoFHandler<dim>::cell_iterator
    cell = dof_handler.begin(),
    endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    { 
				       // Remember the level of the
				       // current cell.
      const unsigned int level = cell->level();
				       // Compute the values specified
				       // by update flags above.
      fe_values.reinit (cell);
      local_matrix = 0;
      
      for (unsigned int q=0; q<n_q_points; ++q)
      {
	  for (unsigned int k=0; k<dofs_per_cell; ++k)
	    {
	      phi_grads_u[k] = fe_values[velocities].gradient (k, q);
	      div_phi_u[k]   = fe_values[velocities].divergence (k, q);
	      phi_p[k]       = fe_values[pressure].value (k, q);
	    }

	  for (unsigned int i=0; i<dofs_per_cell; ++i)
	      for (unsigned int j=0; j<dofs_per_cell; ++j)
		  local_matrix(i,j) += (
                      scalar_product(phi_grads_u[i], phi_grads_u[j])
					- div_phi_u[i] * phi_p[j]
					- phi_p[i] * div_phi_u[j]
					+ phi_p[i] * phi_p[j])
				       * fe_values.JxW(q);     

      }
				      
          cell->get_mg_dof_indices (local_dof_indices);
          boundary_constraints[level]
            .distribute_local_to_global (local_matrix,
                local_dof_indices,
                mg_matrices[level]);

    for (unsigned int i=0; i<dofs_per_cell; ++i)
      for (unsigned int j=0; j<dofs_per_cell; ++j)
	if( !(interface_dofs[level][local_dof_indices[i]]==true && 
	      interface_dofs[level][local_dof_indices[j]]==false))
	  local_matrix(i,j) = 0;

    boundary_interface_constraints[level]
      .distribute_local_to_global (local_matrix,
				   local_dof_indices,
				   mg_interface_matrices[level]);
    }

  mg_A_preconditioner.resize (triangulation.n_levels());
  for (unsigned int level=0; level<triangulation.n_levels(); ++level)
    {
      mg_A_preconditioner[level]
	= std_cxx1x::shared_ptr<typename InnerPreconditioner<dim>::type>(new typename InnerPreconditioner<dim>::type());
      mg_A_preconditioner[level]
	->initialize (mg_matrices[level].block(0,0),
		      typename InnerPreconditioner<dim>::type::AdditionalData());
    }
}

				
template <int dim>
void StokesProblem<dim>::solve_block () 
{
  SparseMatrix<double> pressure_mass_matrix;
  pressure_mass_matrix.reinit(sparsity_pattern.block(1,1));
  pressure_mass_matrix.copy_from(system_matrix.block(1,1));
  system_matrix.block(1,1) = 0;

  SparseILU<double> pmass_preconditioner;
  pmass_preconditioner.initialize (pressure_mass_matrix, 
      SparseILU<double>::AdditionalData());

  InverseMatrix<SparseMatrix<double>,SparseILU<double> >
    m_inverse (pressure_mass_matrix, pmass_preconditioner);

  BlockSchurPreconditioner<typename InnerPreconditioner<dim>::type,
    SparseILU<double> > 
      preconditioner (system_matrix, m_inverse, *A_preconditioner);

  SolverControl solver_control (system_matrix.m(),
      1e-6*system_rhs.l2_norm());
  GrowingVectorMemory<BlockVector<double> > vector_memory;
  SolverGMRES<BlockVector<double> >::AdditionalData gmres_data;
  gmres_data.max_n_tmp_vectors = 100;

  SolverGMRES<BlockVector<double> > gmres(solver_control, vector_memory,
      gmres_data);

  gmres.solve(system_matrix, solution, system_rhs,
      preconditioner);

  constraints.distribute (solution);

  std::cout << " "
    << solver_control.last_step()
    << " block GMRES iterations";
}


template <typename InnerPreconditioner>
class SchurComplementSmoother
{
  public:
    struct AdditionalData 
    {
	const InnerPreconditioner *A_preconditioner;
    };

    void initialize (const BlockSparseMatrix<double> &system_matrix,
		     const AdditionalData            &data);
    
    void vmult (BlockVector<double> &dst,
		const BlockVector<double> &src) const;

    void Tvmult (BlockVector<double> &dst,
		 const BlockVector<double> &src) const;

    void clear ();
    
  private:
    SmartPointer<const BlockSparseMatrix<double> > system_matrix;
    SmartPointer<const InnerPreconditioner>        A_preconditioner;
};


template <typename InnerPreconditioner>
void
SchurComplementSmoother<InnerPreconditioner>::
initialize (const BlockSparseMatrix<double> &system_matrix,
	    const AdditionalData            &data)
{
  this->system_matrix    = &system_matrix;
  this->A_preconditioner = data.A_preconditioner;
}




template <typename InnerPreconditioner>
void
SchurComplementSmoother<InnerPreconditioner>::
vmult (BlockVector<double> &dst,
       const BlockVector<double> &src) const
{
  const InverseMatrix<SparseMatrix<double>,InnerPreconditioner>
    A_inverse (system_matrix->block(0,0), *A_preconditioner);
  Vector<double> tmp (dst.block(0).size());
  
				  
  {
    Vector<double> schur_rhs (dst.block(1).size());
    A_inverse.vmult (tmp, src.block(0));
    system_matrix->block(1,0).vmult (schur_rhs, tmp);
    schur_rhs -= src.block(1);
  
    SchurComplement<InnerPreconditioner>
      schur_complement (*system_matrix, A_inverse);
    
				     // The usual control structures for
				     // the solver call are created...
    SolverControl solver_control (dst.block(1).size(),
				  1e-6*schur_rhs.l2_norm());
    SolverCG<>    cg (solver_control);
    
				    
				   
    SparseILU<double> preconditioner;
    preconditioner.initialize (system_matrix->block(1,1), 
			       SparseILU<double>::AdditionalData());
  
    InverseMatrix<SparseMatrix<double>,SparseILU<double> >
      m_inverse (system_matrix->block(1,1), preconditioner);
    
				    
    cg.solve (schur_complement, dst.block(1), schur_rhs,
	      m_inverse);
  
// no constraints to be taken care of here
    
    std::cout << "  "
	      << solver_control.last_step()
	      << " CG Schur complement iterations in smoother"
	      << std::flush
	      << std::endl;    
  }
    
				 
  {
    system_matrix->block(0,1).vmult (tmp, dst.block(1));
    tmp *= -1;
    tmp += src.block(0);
  
    A_inverse.vmult (dst.block(0), tmp);

// no constraints here either
  }
}



template <typename InnerPreconditioner>
void
SchurComplementSmoother<InnerPreconditioner>::clear ()
{}



template <typename InnerPreconditioner>
void
SchurComplementSmoother<InnerPreconditioner>::
Tvmult (BlockVector<double> &,
	const BlockVector<double> &) const
{
  Assert (false, ExcNotImplemented());
}



template <int dim>
void StokesProblem<dim>::solve_mg () 
{
  assemble_multigrid ();
  typedef PreconditionMG<dim, BlockVector<double>, MGTransferPrebuilt<BlockVector<double> > >
    MGPREC;

  GrowingVectorMemory<BlockVector<double> >  mg_vector_memory;

  MGTransferPrebuilt<BlockVector<double> > mg_transfer(constraints);
  mg_transfer.build_matrices(dof_handler);

  FullMatrix<float> mg_coarse_matrix;
  mg_coarse_matrix.copy_from (mg_matrices[0]);
  MGCoarseGridHouseholder<float, BlockVector<double> > mg_coarse;
  mg_coarse.initialize(mg_coarse_matrix);

  MGMatrix<BlockSparseMatrix<double>, BlockVector<double> >
    mg_matrix(&mg_matrices);

  typedef
    SchurComplementSmoother<typename InnerPreconditioner<dim>::type>
    Smoother;
    
  MGSmootherPrecondition<BlockSparseMatrix<double>,
    Smoother,
    BlockVector<double> >
  mg_smoother(mg_vector_memory);

MGLevelObject<typename Smoother::AdditionalData>
smoother_data (0, triangulation.n_levels());

for (unsigned int level=0; level<triangulation.n_levels(); ++level)
  smoother_data[level].A_preconditioner = mg_A_preconditioner[level].get();

  mg_smoother.initialize(mg_matrices, smoother_data);
  mg_smoother.set_steps(2);

  Multigrid<BlockVector<double> > mg(dof_handler,
      mg_matrix,
      mg_coarse,
      mg_transfer,
      mg_smoother,
      mg_smoother);
  mg.set_debug(3);

  MGPREC  preconditioner(dof_handler, mg, mg_transfer);

  SolverControl solver_control (system_matrix.m(),
      1e-6*system_rhs.l2_norm());
  GrowingVectorMemory<BlockVector<double> > vector_memory;
  SolverGMRES<BlockVector<double> >::AdditionalData gmres_data;
  gmres_data.max_n_tmp_vectors = 100;

  SolverGMRES<BlockVector<double> > gmres(solver_control, vector_memory,
      gmres_data);

//  PreconditionIdentity precondition_identity;
  gmres.solve(system_matrix, solution, system_rhs,
      preconditioner);
  //gmres.solve(system_matrix, solution, system_rhs,
  //    precondition_identity);

  constraints.distribute (solution);

  std::cout << std::endl << " "
    << solver_control.last_step()
    << " block GMRES iterations";
}


template <int dim>
void StokesProblem<dim>::solve_schur () 
{
  const InverseMatrix<SparseMatrix<double>,
                      typename InnerPreconditioner<dim>::type>
    A_inverse (system_matrix.block(0,0), *A_preconditioner);
  Vector<double> tmp (solution.block(0).size());
  
				  
  {
    Vector<double> schur_rhs (solution.block(1).size());
    A_inverse.vmult (tmp, system_rhs.block(0));
    system_matrix.block(1,0).vmult (schur_rhs, tmp);
    schur_rhs -= system_rhs.block(1);
  
    SchurComplement<typename InnerPreconditioner<dim>::type>
      schur_complement (system_matrix, A_inverse);
    
				     // The usual control structures for
				     // the solver call are created...
    SolverControl solver_control (solution.block(1).size(),
				  1e-6*schur_rhs.l2_norm());
    SolverCG<>    cg (solver_control);
    
				    
				   
    SparseILU<double> preconditioner;
    preconditioner.initialize (system_matrix.block(1,1), 
      SparseILU<double>::AdditionalData());
  
    InverseMatrix<SparseMatrix<double>,SparseILU<double> >
      m_inverse (system_matrix.block(1,1), preconditioner);
    
				    
    cg.solve (schur_complement, solution.block(1), schur_rhs,
	      m_inverse);
  
				   
    constraints.distribute (solution);
  
    std::cout << "  "
	      << solver_control.last_step()
	      << " outer CG Schur complement iterations for pressure"
	      << std::flush
	      << std::endl;    
  }
    
				 
  {
    system_matrix.block(0,1).vmult (tmp, solution.block(1));
    tmp *= -1;
    tmp += system_rhs.block(0);
  
    A_inverse.vmult (solution.block(0), tmp);

    constraints.distribute (solution);
  }


}


template <int dim>
void StokesProblem<dim>::find_dofs_on_lower_level (std::vector<std::vector<bool> > &dof_lower_level,
    std::vector<std::vector<bool> > &boundary_dof_lower_level)
{
  const unsigned int   dofs_per_cell   = fe.dofs_per_cell;
  const unsigned int   dofs_per_face   = fe.dofs_per_face;

  std::vector<unsigned int> local_dof_indices (dofs_per_cell);
  std::vector<unsigned int> face_dof_indices (dofs_per_face);

  typename MGDoFHandler<dim>::cell_iterator cell = dof_handler.begin(),
           endc = dof_handler.end();

  for (; cell!=endc; ++cell)
  {
    std::vector<bool> cell_dofs(dofs_per_cell);
    std::vector<bool> boundary_cell_dofs(dofs_per_cell);
    const unsigned int level = cell->level();
    cell->get_mg_dof_indices (local_dof_indices);
    for (unsigned int face_nr=0; 
        face_nr<GeometryInfo<dim>::faces_per_cell; ++face_nr)
    {
      typename DoFHandler<dim>::face_iterator face = cell->face(face_nr);
      if(!cell->at_boundary(face_nr))
      {
        //interior face
        typename MGDoFHandler<dim>::cell_iterator neighbor 
          = cell->neighbor(face_nr);
        // Do refinement face
        // from the coarse side
        if (neighbor->level() < cell->level())
        {
          for(unsigned int j=0; j<dofs_per_face; ++j)
          {
            cell_dofs[fe.face_to_cell_index(j,face_nr)] = true;
          }
        }
      }
      else
      {
        //boundary face
        for(unsigned int inner_face_nr = 1; inner_face_nr<GeometryInfo<dim>::faces_per_cell; ++inner_face_nr)
        {
          const unsigned int neighbor_face_nr = (face_nr+inner_face_nr)%GeometryInfo<dim>::faces_per_cell;
          if(!cell->at_boundary(neighbor_face_nr))
          {
            //other face is interior
            typename MGDoFHandler<dim>::cell_iterator neighbor 
              = cell->neighbor(neighbor_face_nr);
            // Do refinement face
            // from the coarse side
            if (neighbor->level() < cell->level())
            {
              for(unsigned int j=0; j<dofs_per_face; ++j)
              {
                boundary_cell_dofs[fe.face_to_cell_index(j,face_nr)] = true;
              }
            }
          }
        }
      }
    }//faces
    for(unsigned int i=0; i<dofs_per_cell; ++i)
    {
      if(cell_dofs[i])
        dof_lower_level[level][local_dof_indices[i]] = true;

      if(boundary_cell_dofs[i])
        boundary_dof_lower_level[level][local_dof_indices[i]] = true;
    }
  }
}
			
template <int dim>
void
StokesProblem<dim>::output_results (const unsigned int refinement_cycle)  const
{
  std::vector<std::string> solution_names (dim, "velocity");
  solution_names.push_back ("pressure");
  
  std::vector<DataComponentInterpretation::DataComponentInterpretation>
    data_component_interpretation
    (dim, DataComponentInterpretation::component_is_part_of_vector);
  data_component_interpretation
    .push_back (DataComponentInterpretation::component_is_scalar);
      
  DataOut<dim> data_out;
  data_out.attach_dof_handler (dof_handler);  
  data_out.add_data_vector (solution, solution_names,
			    DataOut<dim>::type_dof_data,
			    data_component_interpretation);
  data_out.build_patches ();
  
  std::ostringstream filename;
  filename << "solution-"
           << Utilities::int_to_string (refinement_cycle, 2)
           << ".vtk";

  std::ofstream output (filename.str().c_str());
  data_out.write_vtk (output);
}


			
template <int dim>
void
StokesProblem<dim>::refine_mesh () 
{
  Vector<float> estimated_error_per_cell (triangulation.n_active_cells());

  std::vector<bool> component_mask (dim+1, false);
  component_mask[dim] = true;
  KellyErrorEstimator<dim>::estimate (static_cast<const DoFHandler<dim>&>(dof_handler),
                                      QGauss<dim-1>(degree+1),
                                      typename FunctionMap<dim>::type(),
                                      solution,
                                      estimated_error_per_cell,
                                      component_mask);

  GridRefinement::refine_and_coarsen_fixed_number (triangulation,
                                                   estimated_error_per_cell,
                                                   0.3, 0.0);
  triangulation.execute_coarsening_and_refinement ();
}


				
template <int dim>
void StokesProblem<dim>::run () 
{
  {
    std::vector<unsigned int> subdivisions (dim, 1);
    subdivisions[0] = 4;

    const Point<dim> bottom_left = (dim == 2 ?
				    Point<dim>(0,0) :
				    Point<dim>(0,0,0));
    const Point<dim> top_right   = (dim == 2 ?
				    Point<dim>(1,1) :
				    Point<dim>(1,1,1));
    
    GridGenerator::subdivided_hyper_rectangle (triangulation,
					       subdivisions,
					       bottom_left,
					       top_right);
  }
  
				 
  for (typename Triangulation<dim>::active_cell_iterator
	 cell = triangulation.begin_active();
       cell != triangulation.end(); ++cell)
    for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
      if (cell->face(f)->center()[0] == 1)
	cell->face(f)->set_all_boundary_indicators(1);
  
  
				
  triangulation.refine_global (4-dim);

				  
  for (unsigned int refinement_cycle = 0; refinement_cycle<6;
       ++refinement_cycle)
    {
      std::cout << "Refinement cycle " << refinement_cycle << std::endl;
      
      if (refinement_cycle > 0)
        refine_mesh ();
      
      setup_dofs ();

      std::cout << "   Assembling..." << std::endl << std::flush;
      assemble_system ();      

      std::cout << "   Solving..." << std::flush;
      solve_mg ();
      
      output_results (refinement_cycle);

      std::cout << std::endl;
    }
}


			
int main () 
{
  try
    {
      deallog.depth_console (0);

      StokesProblem<2> flow_problem(1);
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
