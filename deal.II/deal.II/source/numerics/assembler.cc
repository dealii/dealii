/* $Id$ */

#include <numerics/assembler.h>
#include <numerics/base.h>
#include <grid/tria_iterator.h>
#include <grid/tria_iterator.templates.h>
#include <lac/dfmatrix.h>
#include <lac/dvector.h>



template <int dim>
Equation<dim>::Equation (const unsigned int n_equations) :
		n_eq(n_equations) {};



template <int dim>
void Equation<dim>::assemble (dFMatrix          &,
			      dVector           &,
			      const FEValues<dim> &,
			      const typename Triangulation<dim>::cell_iterator &) const {
  Assert (false, ExcPureVirtualFunctionCalled());
};



template <int dim>
void Equation<dim>::assemble (dFMatrix          &,
			      const FEValues<dim> &,
			      const typename Triangulation<dim>::cell_iterator &) const {
  Assert (false, ExcPureVirtualFunctionCalled());
};



template <int dim>
void Equation<dim>::assemble (dVector           &,
			      const FEValues<dim> &,
			      const typename Triangulation<dim>::cell_iterator &) const {
  Assert (false, ExcPureVirtualFunctionCalled());
};




template <int dim>
AssemblerData<dim>::AssemblerData (const DoFHandler<dim>   &dof,
				   const bool               assemble_matrix,
				   const bool               assemble_rhs,
				   ProblemBase<dim>         &problem,
				   const Quadrature<dim>    &quadrature,
				   const FiniteElement<dim> &fe,
				   const FEValues<dim>::UpdateStruct &update_flags) :
		dof(dof),
		assemble_matrix(assemble_matrix),
		assemble_rhs(assemble_rhs),
		problem(problem),
		quadrature(quadrature),
		fe(fe),
		update_flags(update_flags) {};




template <int dim>
Assembler<dim>::Assembler (Triangulation<dim> *tria,
			   const int           level,
			   const int           index,
			   const void         *local_data) :
		DoFCellAccessor<dim> (tria,level,index,
				      &((AssemblerData<dim>*)local_data)->dof),
		cell_matrix (dof_handler->get_selected_fe().total_dofs),
		cell_vector (dVector(dof_handler->get_selected_fe().total_dofs)),
		assemble_matrix (((AssemblerData<dim>*)local_data)->assemble_matrix),
		assemble_rhs (((AssemblerData<dim>*)local_data)->assemble_rhs),
		problem (((AssemblerData<dim>*)local_data)->problem),
		fe(((AssemblerData<dim>*)local_data)->fe),
		fe_values (((AssemblerData<dim>*)local_data)->fe,
			   ((AssemblerData<dim>*)local_data)->quadrature,
			   ((AssemblerData<dim>*)local_data)->update_flags)
{
  Assert ((unsigned int)problem.system_matrix.m() == dof_handler->n_dofs(),
	  ExcInvalidData());
  Assert ((unsigned int)problem.system_matrix.n() == dof_handler->n_dofs(),
	  ExcInvalidData());
  Assert (((AssemblerData<dim>*)local_data)->fe == dof_handler->get_selected_fe(),
	  ExcInvalidData());
  Assert ((unsigned int)problem.right_hand_side.n() == dof_handler->n_dofs(),
	  ExcInvalidData());
};



template <int dim>
void Assembler<dim>::assemble (const Equation<dim> &equation) {
				   // re-init fe values for this cell
  fe_values.reinit (Triangulation<dim>::cell_iterator (tria,
						       present_level,
						       present_index),
		    fe);
  const unsigned int n_dofs = dof_handler->get_selected_fe().total_dofs;

				   // clear cell matrix
  if (assemble_matrix)
    for (unsigned int i=0; i<n_dofs; ++i)
      for (unsigned int j=0; j<n_dofs; ++j)
	cell_matrix(i,j) = 0;
  

				   // clear cell vector
  if (assemble_rhs)
    cell_vector.clear ();
  

				   // fill cell matrix and vector if required
  if (assemble_matrix && assemble_rhs) 
    equation.assemble (cell_matrix, cell_vector, fe_values,
		       Triangulation<dim>::cell_iterator(tria,
							 present_level,
							 present_index));
  else
    if (assemble_matrix)
      equation.assemble (cell_matrix, fe_values,
			 Triangulation<dim>::cell_iterator(tria,
							   present_level,
							   present_index));
    else
      if (assemble_rhs)
	equation.assemble (cell_vector, fe_values,
			   Triangulation<dim>::cell_iterator(tria,
							     present_level,
							     present_index));
      else
	Assert (false, ExcNoAssemblingRequired());


				   // get indices of dofs
  vector<int> dofs;
  get_dof_indices (dofs);
  
				   // distribute cell matrix
  if (assemble_matrix)
    for (unsigned int i=0; i<n_dofs; ++i)
      for (unsigned int j=0; j<n_dofs; ++j)
	problem.system_matrix.add(dofs[i], dofs[j], cell_matrix(i,j));

				   // distribute cell vector
  if (assemble_rhs)
    for (unsigned int j=0; j<n_dofs; ++j)
      problem.right_hand_side(dofs[j]) += cell_vector(j);
};

		



// explicit instantiations
template class Equation<1>;
template class Equation<2>;

template class Assembler<1>;
template class Assembler<2>;


template class AssemblerData<1>;
template class AssemblerData<2>;


template class TriaRawIterator<1,Assembler<1> >;
template class TriaIterator<1,Assembler<1> >;
template class TriaActiveIterator<1,Assembler<1> >;
template class TriaRawIterator<2,Assembler<2> >;
template class TriaIterator<2,Assembler<2> >;
template class TriaActiveIterator<2,Assembler<2> >;
