/* $Id$ */

#include <numerics/assembler.h>
#include <numerics/base.h>
#include <grid/dof_constraints.h>
#include <grid/tria_iterator.h>
#include <basic/data_io.h>
#include <basic/function.h>

#include "../../../mia/control.h"
#include "../../../mia/vectormemory.h"
#include "../../../mia/cg.h"







template <int dim>
ProblemBase<dim>::ProblemBase (Triangulation<dim> *tria,
			       DoFHandler<dim>    *dof) :
		tria(tria),
		dof_handler(dof),
		system_sparsity(1,1,1),        // dummy initialisation, is later reinit'd
		system_matrix()               // dummy initialisation, is later reinit'd
{
  Assert (tria == &dof->get_tria(), ExcDofAndTriaDontMatch());
};



template <int dim>
ProblemBase<dim>::~ProblemBase () {};



template <int dim>
void ProblemBase<dim>::assemble (const Equation<dim>   &equation,
				 const Quadrature<dim> &quadrature,
				 const FiniteElement<dim> &fe) {
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
  solution.reinit (dof_handler->n_dofs());

				   // create assembler
  AssemblerData<dim> data (*dof_handler,
			   true, true, //assemble matrix and rhs
			   *this,
			   quadrature,
			   fe);
  TriaActiveIterator<dim, Assembler<dim> > assembler (tria,
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
};




template <int dim>
void ProblemBase<dim>::solve () {
  int    max_iter  = 4000;
  double tolerance = 1.e-16;
  
  Control                          control1(max_iter,tolerance);
  PrimitiveVectorMemory<dVector>   memory(right_hand_side.n());
  CG<dSMatrix,dVector>             cg(control1,memory);

				   // solve
  cg (system_matrix, solution, right_hand_side);
				   // distribute solution
  constraints.distribute (solution);
};



/*
template <int dim>
void ProblemBase<dim>::integrate_L1_difference (const Function<dim>      &exact_solution,
						vector<double>           &difference,
						const Quadrature<dim>    &q,
						const FiniteElement<dim> &fe) const {
  Assert (fe == dof_handler->get_selected_fe(), ExcInvalidFE());

  difference.erase (difference.begin(), difference.end());
  difference.reserve (tria->n_cells());

				   // loop over all cells
  DoFHandler<dim>::active_dof_iterator cell = dof_handler->begin_active(),
				       endc = dof_handler->end();
  for (; cell != end; ++cell) 
    {
      double diff=0;
      
      difference.push_back (diff);
    };
};
*/



template <int dim>
void ProblemBase<dim>::fill_data (DataOut<dim> &out) const {
  out.clear_data_vectors ();
  out.attach_dof_handler (*dof_handler);

  pair<char*,char*> solution_name = get_solution_name ();
  out.add_data_vector (solution,
		       solution_name.first, solution_name.second);
};


/*
template <int dim>
pair<char*,char*> ProblemBase<dim>::get_solution_name () const {
  return make_pair("solution", "");
};
*/

  
  


// explicit instantiations
template class ProblemBase<1>;
template class ProblemBase<2>;
