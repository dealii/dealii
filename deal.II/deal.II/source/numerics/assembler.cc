/* $Id$ */

#include <numerics/assembler.h>
#include <numerics/base.h>
#include <grid/tria_iterator.h>
#include <lac/dfmatrix.h>
#include <lac/dvector.h>


extern TriaIterator<1,CellAccessor<1> > __dummy127; // do this to calm down gcc2.7,
extern TriaIterator<2,CellAccessor<2> > __dummy128; // wait for gcc2.8


template <int dim>
Equation<dim>::Equation (const unsigned int n_equations) :
		n_eq(n_equations) {};



void Equation<1>::assemble (dFMatrix          &,
			    vector<dVector>   &,
			    const FEValues<1> &,
			    const Triangulation<1>::cell_iterator &) const {
  Assert (false, ExcPureVirtualFunctionCalled());
};



void Equation<2>::assemble (dFMatrix          &,
			    vector<dVector>   &,
			    const FEValues<2> &,
			    const Triangulation<2>::cell_iterator &) const {
  Assert (false, ExcPureVirtualFunctionCalled());
};



void Equation<1>::assemble (dFMatrix          &,
			    const FEValues<1> &,
			    const Triangulation<1>::cell_iterator &) const {
  Assert (false, ExcPureVirtualFunctionCalled());
};



void Equation<2>::assemble (dFMatrix          &,
			    const FEValues<2> &,
			    const Triangulation<2>::cell_iterator &) const {
  Assert (false, ExcPureVirtualFunctionCalled());
};



void Equation<1>::assemble (vector<dVector>   &,
			    const FEValues<1> &,
			    const Triangulation<1>::cell_iterator &) const {
  Assert (false, ExcPureVirtualFunctionCalled());
};



void Equation<2>::assemble (vector<dVector>   &,
			    const FEValues<2> &,
			    const Triangulation<2>::cell_iterator &) const {
  Assert (false, ExcPureVirtualFunctionCalled());
};



/*
template <int dim>
Assembler<dim>::AssemblerData::AssemblerData (DoFHandler<dim> *dof,
					      const bool assemble_matrix,
					      const bool assemble_rhs,
					      const unsigned int n_rhs,
					      ProblemBase<dim> *problem,
					      const Quadrature<dim> *quadrature) :
		dof(dof), assemble_matrix(assemble_matrix),
		assemble_rhs(assemble_rhs), n_rhs(n_rhs),
		problem(problem), quadrature(quadrature) {};
		*/



template <int dim>
Assembler<dim>::Assembler (Triangulation<dim> *tria,
			   const int           level,
			   const int           index,
			   const void         *local_data) :
		DoFCellAccessor<dim> (tria,level,index,
				      ((AssemblerData*)local_data)->dof),
		cell_matrix (dof_handler->get_selected_fe().total_dofs),
		cell_vectors (((AssemblerData*)local_data)->n_rhs,
			      dVector(dof_handler->get_selected_fe().total_dofs)),
		assemble_matrix (((AssemblerData*)local_data)->assemble_matrix),
		assemble_rhs (((AssemblerData*)local_data)->assemble_rhs),
		problem (((AssemblerData*)local_data)->problem),
		fe_values (dof_handler->get_selected_fe(),
			   *((AssemblerData*)local_data)->quadrature)
{
  Assert (((AssemblerData*)local_data)->dof != 0,
	  ExcInvalidData());
  Assert (((AssemblerData*)local_data)->n_rhs != 0,
	  ExcInvalidData());
  Assert (((AssemblerData*)local_data)->problem != 0,
	  ExcInvalidData());

  Assert ((unsigned int)problem->system_matrix.m() == dof_handler->n_dofs(),
	  ExcInvalidData());
  Assert ((unsigned int)problem->system_matrix.n() == dof_handler->n_dofs(),
	  ExcInvalidData());
  Assert (problem->right_hand_sides.size() == cell_vectors.size(),
	  ExcInvalidData());
  for (unsigned int i=0; i<cell_vectors.size(); ++i)
    Assert ((unsigned int)problem->right_hand_sides[i]->n() == dof_handler->n_dofs(),
	    ExcInvalidData());
};



template <int dim>
void Assembler<dim>::assemble (const Equation<dim> &equation) {
				   // re-init fe values for this cell
  fe_values.reinit (Triangulation<dim>::cell_iterator (tria,
						       present_level,
						       present_index),
		    dof_handler->get_selected_fe());
  
  if (assemble_matrix)
				     // clear cell matrix
    for (unsigned int i=0; i<dof_handler->get_selected_fe().total_dofs; ++i)
      for (unsigned int j=0; j<dof_handler->get_selected_fe().total_dofs; ++j)
	cell_matrix(i,j) = 0;
  

  if (assemble_rhs)
				     // clear cell vector
    for (unsigned int vec=0; vec<cell_vectors.size(); ++vec)
      for (unsigned int j=0; j<dof_handler->get_selected_fe().total_dofs; ++j)
	cell_vectors[vec](j) = 0;
  

				   // fill cell matrix and vector if required
  if (assemble_matrix && assemble_rhs) 
    equation.assemble (cell_matrix, cell_vectors, fe_values,
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
	equation.assemble (cell_vectors, fe_values,
			   Triangulation<dim>::cell_iterator(tria,
							     present_level,
							     present_index));
      else
	Assert (false, ExcNoAssemblingRequired());


				   // get indices of dofs
  vector<int> dofs;
  dof_indices (dofs);
  
  if (assemble_matrix)
				     // distribute cell matrix
    for (unsigned int i=0; i<dof_handler->get_selected_fe().total_dofs; ++i)
      for (unsigned int j=0; j<dof_handler->get_selected_fe().total_dofs; ++j)
	problem->system_matrix.set(dofs[i], dofs[j], cell_matrix(i,j));

  if (assemble_rhs)
				     // distribute cell vector
    for (unsigned int vec=0; vec<cell_vectors.size(); ++vec)
      for (unsigned int j=0; j<dof_handler->get_selected_fe().total_dofs; ++j)
	(*problem->right_hand_sides[vec])(dofs[j]) = cell_vectors[vec](j);
};

		



// explicit instantiations
template class Equation<1>;
template class Equation<2>;

template class Assembler<1>;
template class Assembler<2>;
