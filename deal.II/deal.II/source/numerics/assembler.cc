/* $Id$ */
/* Copyright W. Bangerth, University of Heidelberg, 1998 */


#include <numerics/assembler.h>
#include <grid/tria_iterator.h>
#include <grid/tria_iterator.templates.h>
#include <fe/fe.h>
#include <lac/fullmatrix.h>
#include <lac/vector.h>
#include <lac/sparsematrix.h>


template <int dim>
Assembler<dim>::AssemblerData::AssemblerData (const DoFHandler<dim>    &dof,
					      const bool                assemble_matrix,
					      const bool                assemble_rhs,
					      SparseMatrix<double>     &matrix,
					      Vector<double>           &rhs_vector,
					      const Quadrature<dim>    &quadrature,
					      const UpdateFlags        &update_flags) :
		dof(dof),
		assemble_matrix(assemble_matrix),
		assemble_rhs(assemble_rhs),
		matrix(matrix),
		rhs_vector(rhs_vector),
		quadrature(quadrature),
		update_flags(update_flags)
{};




template <int dim>
Assembler<dim>::Assembler (Triangulation<dim>  *tria,
			   const int            level,
			   const int            index,
			   const AssemblerData *local_data) :
		DoFCellAccessor<dim> (tria,level,index, &local_data->dof),
		cell_matrix (dof_handler->get_fe().total_dofs),
		cell_vector (Vector<double>(dof_handler->get_fe().total_dofs)),
		assemble_matrix (local_data->assemble_matrix),
		assemble_rhs (local_data->assemble_rhs),
		matrix(local_data->matrix),
		rhs_vector(local_data->rhs_vector),
		fe_values (dof_handler->get_fe(),
			   local_data->quadrature,
			   local_data->update_flags,
			   tria->get_boundary())
{
  Assert (!assemble_matrix || (matrix.m() == dof_handler->n_dofs()),
	  ExcInvalidData());
  Assert (!assemble_matrix || (matrix.n() == dof_handler->n_dofs()),
	  ExcInvalidData());
  Assert (!assemble_rhs || (rhs_vector.size()==dof_handler->n_dofs()),
	  ExcInvalidData());
};



template <int dim>
void Assembler<dim>::assemble (const Equation<dim> &equation) {
				   // re-init fe values for this cell
  fe_values.reinit (DoFHandler<dim>::cell_iterator (*this));
  const unsigned int n_dofs = dof_handler->get_fe().total_dofs;

  if (assemble_matrix)
    cell_matrix.clear ();
  if (assemble_rhs)
    cell_vector.clear ();
  

				   // fill cell matrix and vector if required
  DoFHandler<dim>::cell_iterator this_cell (*this);
  if (assemble_matrix && assemble_rhs) 
    equation.assemble (cell_matrix, cell_vector, fe_values, this_cell);
  else
    if (assemble_matrix)
      equation.assemble (cell_matrix, fe_values, this_cell);
    else
      if (assemble_rhs)
	equation.assemble (cell_vector, fe_values, this_cell);
      else
	Assert (false, ExcNoAssemblingRequired());


				   // get indices of dofs
  vector<int> dofs (n_dofs);
  get_dof_indices (dofs);

				   // one could use the
				   // #distribute_local_to_global# functions
				   // here, but they would require getting the
				   // dof indices twice, so we leave it the
				   // way it was originally programmed.
  
				   // distribute cell matrix
  if (assemble_matrix)
    for (unsigned int i=0; i<n_dofs; ++i)
      for (unsigned int j=0; j<n_dofs; ++j)
	matrix.add(dofs[i], dofs[j], cell_matrix(i,j));

				   // distribute cell vector
  if (assemble_rhs)
    for (unsigned int j=0; j<n_dofs; ++j)
      rhs_vector(dofs[j]) += cell_vector(j);
};

		



// explicit instantiations
template class Equation<deal_II_dimension>;
template class Assembler<deal_II_dimension>;

template class TriaRawIterator<deal_II_dimension,Assembler<deal_II_dimension> >;
template class TriaIterator<deal_II_dimension,Assembler<deal_II_dimension> >;
template class TriaActiveIterator<deal_II_dimension,Assembler<deal_II_dimension> >;
