/* $Id$ */

#include <numerics/assembler.h>
#include <numerics/base.h>
#include <grid/dof_constraints.h>
#include <grid/tria_iterator.h>


extern TriaActiveIterator<1,CellAccessor<1> > __dummy1233; // for gcc2.8
extern TriaActiveIterator<2,CellAccessor<2> > __dummy1234;



template <int dim>
ProblemBase<dim>::ProblemBase (Triangulation<dim> *tria,
			       DoFHandler<dim>    *dof,
			       const unsigned int  n_rhs) :
		tria(tria),
		dof_handler(dof),
		system_sparsity(1,1,1),        // dummy initialisation, is later reinit'd
		system_matrix(system_sparsity),// dummy initialisation, is later reinit'd
		right_hand_sides (n_rhs, (dVector*)0)	
{
  Assert (tria == &dof->get_tria(), ExcDofAndTriaDontMatch());

  for (unsigned int i=0; i<n_rhs; ++i) 
    {
      right_hand_sides[i] = new dVector;
      Assert (right_hand_sides[i] != 0, ExcNoMemory());
    };
};




template <int dim>
void ProblemBase<dim>::assemble (const Equation<dim>   &equation,
				 const Quadrature<dim> &quadrature) {
  system_sparsity.reinit (dof_handler->n_dofs(),
			  dof_handler->n_dofs(),
			  dof_handler->max_couplings_between_dofs());
  solution.reinit (dof_handler->n_dofs());
  for (unsigned int i=0; i<right_hand_sides.size(); ++i)
    right_hand_sides[i]->reinit (dof_handler->n_dofs());
  
				   // make up sparsity pattern and
				   // compress with constraints
  ConstraintMatrix constraints;
  
  dof_handler->make_constraint_matrix (constraints);
  dof_handler->make_sparsity_pattern (system_sparsity);
  constraints.condense (system_sparsity);

				   // reinite system matrix
  system_matrix.reinit (system_sparsity);

				   // create assembler
  Assembler<dim>::AssemblerData data (dof_handler,
				      true, true, //assemble matrix and rhs
				      right_hand_sides.size(),   //one rhs
				      this,
				      &quadrature);
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
};






// explicit instantiations
template class ProblemBase<1>;
template class ProblemBase<2>;
