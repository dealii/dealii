/* $Id$ */

#include <numerics/assembler.h>
#include <numerics/base.h>



template <int dim>
ProblemBase<dim>::ProblemBase (Triangulation<dim> *tria,
			       DoFHandler<dim>    *dof) :
		tria(tria), dof_handler(dof) {};




template <int dim>
void ProblemBase<dim>::assemble (const Equation<dim>   &equation,
				 const Quadrature<dim> &quadrature) {
  ConstraintMatrix constraints;

  system_sparsity.reinit (dof_handler->n_dofs(),
			  dof_handler->n_dofs(),
			  dof_handler->max_coupling_between_dofs());
  
  dof_handler->make_constraint_matrix (constraints);
  dof_handler->make_sparsity_pattern (system_sparsity);
  constraints.condense (system_sparsity);

  system_matrix.reinit (system_sparsity);

  Assembler<dim>::AssemblerData data = {dof_handler,
					true, true,
					1,
					this,
					&quadrature};
  
  TriaActiveIterator<Assembler<dim> > assembler (tria,
						 tria->begin_active()->level(),
						 tria->begin_active()->index(),
						 &data);

  do 
    {
      assembler->assemble (equation);
    }
  while ((++assembler).state() == valid);
};
