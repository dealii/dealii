/* $Id$ */

#include <numerics/assembler.h>
#include <numerics/base.h>
#include <grid/dof_constraints.h>
#include <grid/tria_iterator.h>
#include <basic/data_io.h>


#include "../../../mia/vectormemory.h"
#include "../../../mia/control.h"
#include "../../../mia/cg.h"

extern TriaActiveIterator<1,CellAccessor<1> > __dummy1233; // for gcc2.7
extern TriaActiveIterator<2,CellAccessor<2> > __dummy1234;



template <int dim>
ProblemBase<dim>::ProblemBase (Triangulation<dim> *tria,
			       DoFHandler<dim>    *dof,
			       const unsigned int  n_rhs) :
		tria(tria),
		dof_handler(dof),
		system_sparsity(1,1,1),        // dummy initialisation, is later reinit'd
		system_matrix(),               // dummy initialisation, is later reinit'd
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
				 const Quadrature<dim> &quadrature,
				 const FiniteElement<dim> &fe) {
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
  AssemblerData<dim> data (*dof_handler,
			   true, true, //assemble matrix and rhs
			   right_hand_sides.size(),   //one rhs
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
};




template <int dim>
void ProblemBase<dim>::solve () {
  int    max_iter  = 4000;
  double tolerance = 1.e-14;
  
  Control                          control1(max_iter,tolerance);
  PrimitiveVectorMemory<dVector>   memory(right_hand_sides[0]->n());
  CG<dSMatrix,dVector>             cg(control1,memory);

  cg (system_matrix, solution, *right_hand_sides[0]);
};



template <int dim>
void ProblemBase<dim>::fill_data (DataOut<dim> &out) const {
  out.clear_data_vectors ();
  out.attach_dof_handler (*dof_handler);
  out.add_data_vector (solution, "solution");
};

  
  


// explicit instantiations
template class ProblemBase<1>;
template class ProblemBase<2>;
