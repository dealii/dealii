/* $Id$ */
/* Copyright W. Bangerth, University of Heidelberg, 1998 */


#include <numerics/assembler.h>
#include <numerics/base.h>
#include <numerics/matrices.h>
#include <numerics/vectors.h>
#include <grid/dof_constraints.h>
#include <grid/tria_iterator.h>
#include <basic/data_io.h>
#include <basic/function.h>
#include <fe/fe.h>
#include <fe/quadrature.h>
#include <lac/dvector.h>
#include <lac/solver_cg.h>
#include <lac/vector_memory.h>

#include <map>
#include <numeric>
#include <algorithm>
#include <cmath>




template <int dim>
ProblemBase<dim>::ProblemBase () :
		tria(0),
		dof_handler(0),
		system_sparsity(),        // dummy initialisation, is later reinit'd
		system_matrix()           // dummy initialisation, is later reinit'd
{};



template <int dim>
void ProblemBase<dim>::set_tria_and_dof (Triangulation<dim> *t,
					 DoFHandler<dim>    *d) {
  tria        = t;
  dof_handler = d;

  Assert ((tria!=0) && (dof_handler!=0), ExcNoTriaSelected());
  Assert (tria == &dof_handler->get_tria(), ExcDofAndTriaDontMatch());
};



template <int dim>
void ProblemBase<dim>::clear () {
  if (tria)        { delete tria;         tria        = 0; };
  if (dof_handler) { delete dof_handler;  dof_handler = 0; };
  system_sparsity.reinit (0,0,1);
  system_matrix.clear ();
  right_hand_side.reinit (1);
  solution.reinit (1);
  constraints.clear ();
};



template <int dim>
ProblemBase<dim>::~ProblemBase () {};



template <int dim>
void ProblemBase<dim>::assemble (const Equation<dim>      &equation,
				 const Quadrature<dim>    &quadrature,
				 const FiniteElement<dim> &fe,
				 const UpdateFlags         update_flags,
				 const FunctionMap        &dirichlet_bc,
				 const Boundary<dim>      &boundary) {
  Assert ((tria!=0) && (dof_handler!=0), ExcNoTriaSelected());
  
  system_sparsity.reinit (dof_handler->n_dofs(),
			  dof_handler->n_dofs(),
			  dof_handler->max_couplings_between_dofs());
  right_hand_side.reinit (dof_handler->n_dofs());
  
				   // make up sparsity pattern and
				   // compress with constraints
  constraints.clear ();
  dof_handler->make_constraint_matrix (constraints);
  dof_handler->make_sparsity_pattern (system_sparsity);
  constraints.condense (system_sparsity);

				   // reinite system matrix
  system_matrix.reinit (system_sparsity);
				   // reinit solution vector, preset
				   // with zeroes.
  solution.reinit (dof_handler->n_dofs());
  
				   // create assembler
  AssemblerData<dim> data (*dof_handler,
			   true, true, //assemble matrix and rhs
			   system_matrix,
			   right_hand_side,
			   quadrature,
			   fe,
			   update_flags,
			   boundary);
  active_assemble_iterator assembler (tria,
				      tria->begin_active()->level(),
				      tria->begin_active()->index(),
				      &data);
				   // loop over all cells, fill matrix and rhs
  do 
    {
      assembler->assemble (equation);
    }
  while ((++assembler).state() == valid);

				   // condense system matrix in-place
  constraints.condense (system_matrix);

				   // condense right hand side in-place
  constraints.condense (right_hand_side);

				   // apply Dirichlet bc as described
				   // in the docs
  map<int, double> boundary_value_list;
  VectorTools<dim>::interpolate_boundary_values (*dof_handler,
						 dirichlet_bc, fe, boundary,
						 boundary_value_list);
  MatrixTools<dim>::apply_boundary_values (boundary_value_list,
					   system_matrix, solution,
					   right_hand_side);  
};




template <int dim>
void ProblemBase<dim>::solve () {
  Assert ((tria!=0) && (dof_handler!=0), ExcNoTriaSelected());
  
  SolverControl                    control(4000, 1e-16);
  PrimitiveVectorMemory<dVector>   memory;
  SolverCG<dSMatrix,dVector>       cg(control,memory);

				   // solve
  cg.solve (system_matrix, solution, right_hand_side);
				   // distribute solution
  constraints.distribute (solution);
};




template <int dim>
void ProblemBase<dim>::fill_data (DataOut<dim> &out) const {
  Assert ((tria!=0) && (dof_handler!=0), ExcNoTriaSelected());
  
  out.clear_data_vectors ();
  out.attach_dof_handler (*dof_handler);

  pair<char*,char*> solution_name = get_solution_name ();
  out.add_data_vector (solution,
		       solution_name.first, solution_name.second);
};



template <int dim>
pair<char*,char*> ProblemBase<dim>::get_solution_name () const {
  return pair<char*,char*>("solution", "<dimensionless>");
};




  


// explicit instantiations
template class ProblemBase<deal_II_dimension>;
