/* $Id$ */
/* Copyright W. Bangerth, University of Heidelberg, 1998 */



#include <grid/dof_accessor.h>
#include <grid/dof_accessor.templates.h>
#include <grid/dof.h>
#include <grid/dof_levels.h>
#include <grid/tria_iterator.h>
#include <grid/tria_iterator.templates.h>
#include <fe/fe.h>

#include <lac/vector.h>
#include <lac/fullmatrix.h>
#include <lac/sparsematrix.h>

#include <vector>





/*------------------------- Functions: DoFLineAccessor -----------------------*/



template <int dim>
void DoFObjectAccessor<1, dim>::set_dof_index (const unsigned int i,
					       const int index) const {
  Assert (dof_handler != 0, DoFAccessor<dim>::ExcInvalidObject());
				   // make sure a FE has been selected
				   // and enough room was reserved
  Assert (dof_handler->selected_fe != 0, DoFAccessor<dim>::ExcInvalidObject());
  Assert (i<dof_handler->selected_fe->dofs_per_line,
	  ExcIndexRange (i, 0, dof_handler->selected_fe->dofs_per_line));

  dof_handler->levels[present_level]
    ->line_dofs[present_index*dof_handler->selected_fe->dofs_per_line+i] = index;
};




template <int dim>
void DoFObjectAccessor<1, dim>::set_vertex_dof_index (const unsigned int vertex,
						      const unsigned int i,
						      const int index) const {
  Assert (dof_handler != 0, ExcInvalidObject());
  Assert (dof_handler->selected_fe != 0, ExcInvalidObject());
  Assert (vertex<2, ExcIndexRange (i,0,2));
  Assert (i<dof_handler->selected_fe->dofs_per_vertex,
	  ExcIndexRange (i, 0, dof_handler->selected_fe->dofs_per_vertex));

  const unsigned int dof_number = (vertex_index(vertex) *
				   dof_handler->selected_fe->dofs_per_vertex +
				   i);
  dof_handler->vertex_dofs[dof_number] = index;
};



template <int dim>
void DoFObjectAccessor<1, dim>::
distribute_local_to_global (const Vector<double> &local_source,
			    Vector<double>       &global_destination) const {
  Assert (dof_handler != 0, ExcInvalidObject());
  Assert (dof_handler->selected_fe != 0, ExcInvalidObject());
  Assert (local_source.size() == (2*dof_handler->get_fe().dofs_per_vertex +
				  dof_handler->get_fe().dofs_per_line),
	  ExcVectorDoesNotMatch());
  Assert (dof_handler->n_dofs() == global_destination.size(),
	  ExcVectorDoesNotMatch());

  const unsigned int n_dofs = local_source.size();
  
				   // get indices of dofs
  vector<int> dofs (n_dofs);
  get_dof_indices (dofs);
  
				   // distribute cell vector
  for (unsigned int j=0; j<n_dofs; ++j)
    global_destination(dofs[j]) += local_source(j);
};



template <int dim>
void DoFObjectAccessor<1, dim>::
distribute_local_to_global (const FullMatrix<double> &local_source,
			    SparseMatrix<double>     &global_destination) const {
  Assert (dof_handler != 0, ExcInvalidObject());
  Assert (dof_handler->selected_fe != 0, ExcInvalidObject());
  Assert (local_source.m() == (2*dof_handler->get_fe().dofs_per_vertex +
			       dof_handler->get_fe().dofs_per_line),
	  ExcVectorDoesNotMatch());
  Assert (local_source.m() == local_source.n(),
	  ExcMatrixDoesNotMatch());
  Assert (dof_handler->n_dofs() == global_destination.m(),
	  ExcMatrixDoesNotMatch());
  Assert (global_destination.m() == global_destination.n(),
	  ExcMatrixDoesNotMatch());

  const unsigned int n_dofs = local_source.m();

				   // get indices of dofs
  vector<int> dofs (n_dofs);
  get_dof_indices (dofs);
  
				   // distribute cell matrix
  for (unsigned int i=0; i<n_dofs; ++i)
    for (unsigned int j=0; j<n_dofs; ++j)
      global_destination.add(dofs[i], dofs[j], local_source(i,j));
};



/*------------------------- Functions: DoFQuadAccessor -----------------------*/


template <int dim>
void DoFObjectAccessor<2, dim>::set_dof_index (const unsigned int i,
					       const int index) const {
  Assert (dof_handler != 0, ExcInvalidObject());
				   // make sure a FE has been selected
				   // and enough room was reserved
  Assert (dof_handler->selected_fe != 0, ExcInvalidObject());
  Assert (i<dof_handler->selected_fe->dofs_per_quad,
	  ExcIndexRange (i, 0, dof_handler->selected_fe->dofs_per_quad));

  dof_handler->levels[present_level]
    ->quad_dofs[present_index*dof_handler->selected_fe->dofs_per_quad+i] = index;
};



template <int dim>
void DoFObjectAccessor<2, dim>::set_vertex_dof_index (const unsigned int vertex,
							   const unsigned int i,
							   const int index) const {
  Assert (dof_handler != 0, ExcInvalidObject());
  Assert (dof_handler->selected_fe != 0, ExcInvalidObject());
  Assert (vertex<4, ExcIndexRange (i,0,4));
  Assert (i<dof_handler->selected_fe->dofs_per_vertex,
	  ExcIndexRange (i, 0, dof_handler->selected_fe->dofs_per_vertex));

  const unsigned int dof_number = (vertex_index(vertex) *
				   dof_handler->selected_fe->dofs_per_vertex +
				   i);
  dof_handler->vertex_dofs[dof_number] = index;
};



template <int dim>
void DoFObjectAccessor<2, dim>::
distribute_local_to_global (const Vector<double> &local_source,
			    Vector<double>       &global_destination) const {
  Assert (dof_handler != 0, ExcInvalidObject());
  Assert (dof_handler->selected_fe != 0, ExcInvalidObject());
  Assert (local_source.size() == (4*dof_handler->get_fe().dofs_per_vertex +
				  4*dof_handler->get_fe().dofs_per_line +
				  dof_handler->get_fe().dofs_per_quad),
	  ExcVectorDoesNotMatch());
  Assert (dof_handler->n_dofs() == global_destination.size(),
	  ExcVectorDoesNotMatch());

  const unsigned int n_dofs = local_source.size();

				   // get indices of dofs
  vector<int> dofs (n_dofs);
  get_dof_indices (dofs);
  
				   // distribute cell vector
  for (unsigned int j=0; j<n_dofs; ++j)
    global_destination(dofs[j]) += local_source(j);
};



template <int dim>
void DoFObjectAccessor<2, dim>::
distribute_local_to_global (const FullMatrix<double> &local_source,
			    SparseMatrix<double>     &global_destination) const {
  Assert (dof_handler != 0, ExcInvalidObject());
  Assert (dof_handler->selected_fe != 0, ExcInvalidObject());
  Assert (local_source.m() == (4*dof_handler->get_fe().dofs_per_vertex +
			       4*dof_handler->get_fe().dofs_per_line +
			       dof_handler->get_fe().dofs_per_quad),
	  ExcMatrixDoesNotMatch());
  Assert (local_source.m() == local_source.n(),
	  ExcMatrixDoesNotMatch());
  Assert (dof_handler->n_dofs() == global_destination.m(),
	  ExcMatrixDoesNotMatch());
  Assert (global_destination.m() == global_destination.n(),
	  ExcMatrixDoesNotMatch());
  
  const unsigned int n_dofs = local_source.m();

				   // get indices of dofs
  vector<int> dofs (n_dofs);
  get_dof_indices (dofs);
  
				   // distribute cell matrix
  for (unsigned int i=0; i<n_dofs; ++i)
    for (unsigned int j=0; j<n_dofs; ++j)
      global_destination.add(dofs[i], dofs[j], local_source(i,j));
};




/*------------------------- Functions: DoFObjectAccessor -----------------------*/


template <int dim>
void DoFObjectAccessor<3, dim>::set_dof_index (const unsigned int i,
						    const int index) const {
  Assert (dof_handler != 0, ExcInvalidObject());
				   // make sure a FE has been selected
				   // and enough room was reserved
  Assert (dof_handler->selected_fe != 0, ExcInvalidObject());
  Assert (i<dof_handler->selected_fe->dofs_per_hex,
	  ExcIndexRange (i, 0, dof_handler->selected_fe->dofs_per_hex));

  dof_handler->levels[present_level]
    ->hex_dofs[present_index*dof_handler->selected_fe->dofs_per_hex+i] = index;
};



template <int dim>
void DoFObjectAccessor<3, dim>::set_vertex_dof_index (const unsigned int vertex,
							   const unsigned int i,
							   const int index) const {
  Assert (dof_handler != 0, ExcInvalidObject());
  Assert (dof_handler->selected_fe != 0, ExcInvalidObject());
  Assert (vertex<8, ExcIndexRange (i,0,8));
  Assert (i<dof_handler->selected_fe->dofs_per_vertex,
	  ExcIndexRange (i, 0, dof_handler->selected_fe->dofs_per_vertex));

  const unsigned int dof_number = (vertex_index(vertex) *
				   dof_handler->selected_fe->dofs_per_vertex +
				   i);
  dof_handler->vertex_dofs[dof_number] = index;
};



template <int dim>
void DoFObjectAccessor<3, dim>::
distribute_local_to_global (const Vector<double> &local_source,
			    Vector<double>       &global_destination) const {
  Assert (dof_handler != 0, ExcInvalidObject());
  Assert (dof_handler->selected_fe != 0, ExcInvalidObject());
  Assert (local_source.size() == (8*dof_handler->get_fe().dofs_per_vertex +
				  12*dof_handler->get_fe().dofs_per_line +
				  6*dof_handler->get_fe().dofs_per_quad +
				  dof_handler->get_fe().dofs_per_hex),
	  ExcVectorDoesNotMatch());
  Assert (dof_handler->n_dofs() == global_destination.size(),
	  ExcVectorDoesNotMatch());

  const unsigned int n_dofs = local_source.size();

				   // get indices of dofs
  vector<int> dofs (n_dofs);
  get_dof_indices (dofs);
  
				   // distribute cell vector
  for (unsigned int j=0; j<n_dofs; ++j)
    global_destination(dofs[j]) += local_source(j);
};



template <int dim>
void DoFObjectAccessor<3, dim>::
distribute_local_to_global (const FullMatrix<double> &local_source,
			    SparseMatrix<double>     &global_destination) const {
  Assert (dof_handler != 0, ExcInvalidObject());
  Assert (dof_handler->selected_fe != 0, ExcInvalidObject());
  Assert (local_source.m() == (8*dof_handler->get_fe().dofs_per_vertex +
			       12*dof_handler->get_fe().dofs_per_line +
			       6*dof_handler->get_fe().dofs_per_quad +
			       dof_handler->get_fe().dofs_per_hex),
	  ExcMatrixDoesNotMatch());
  Assert (local_source.m() == local_source.n(),
	  ExcMatrixDoesNotMatch());
  Assert (dof_handler->n_dofs() == global_destination.m(),
	  ExcMatrixDoesNotMatch());
  Assert (global_destination.m() == global_destination.n(),
	  ExcMatrixDoesNotMatch());
  
  const unsigned int n_dofs = local_source.m();

				   // get indices of dofs
  vector<int> dofs (n_dofs);
  get_dof_indices (dofs);
  
				   // distribute cell matrix
  for (unsigned int i=0; i<n_dofs; ++i)
    for (unsigned int j=0; j<n_dofs; ++j)
      global_destination.add(dofs[i], dofs[j], local_source(i,j));
};






/*------------------------- Functions: DoFCellAccessor -----------------------*/



#if deal_II_dimension == 1

template <>
TriaIterator<1, DoFObjectAccessor<0,1> >
DoFCellAccessor<1>::face (const unsigned int) const
{
  Assert (false, ExcNotUsefulForThisDimension());
  return TriaIterator<1, DoFObjectAccessor<0,1> >();
};



/* this function looks like a template, but actually isn't one, it
   will only work for the dimension selected by the preprocessor
   guard above.
   
   the correct thing to do would be

   template <>
   template <typename number>
   void
   DoFCellAccessor<1>::get... (vector<number>...)

   but this does not work, at least not at present using egcs1.1.1. we therefore
   separate the different implementations for the different dimensions using
   the preprocecssor and double check using an assertion in the function body.
*/
template <int dim>
template <typename number>
void
DoFCellAccessor<dim>::get_dof_values (const Vector<number> &values,
				      Vector<number>       &local_values) const {
  Assert (dim==1, ExcInternalError());
  
  Assert (dof_handler != 0, DoFAccessor<1>::ExcInvalidObject());
  Assert (&dof_handler->get_fe() != 0, DoFAccessor<1>::ExcInvalidObject());
  Assert (local_values.size() == dof_handler->get_fe().total_dofs,
	  DoFAccessor<1>::ExcVectorDoesNotMatch());
  Assert (values.size() == dof_handler->n_dofs(),
	  DoFAccessor<1>::ExcVectorDoesNotMatch());
  Assert (active(), ExcNotActive());
  
  const unsigned int dofs_per_vertex = dof_handler->get_fe().dofs_per_vertex,
		     dofs_per_line   = dof_handler->get_fe().dofs_per_line;
  typename vector<number>::iterator next_local_value=local_values.begin();
  for (unsigned int vertex=0; vertex<2; ++vertex)
    for (unsigned int d=0; d<dofs_per_vertex; ++d)
      *next_local_value++ = values(vertex_dof_index(vertex,d));
  for (unsigned int d=0; d<dofs_per_line; ++d)
    *next_local_value++ = values(dof_index(d));

  Assert (next_local_value == local_values.end(),
	  ::ExcInternalError());
};



/* this function looks like a template, but actually isn't one, it
   will only work for the dimension selected by the preprocessor
   guard above.
   
   the correct thing to do would be

   template <>
   template <typename number>
   void
   DoFCellAccessor<1>::get... (vector<number>...)

   but this does not work, at least not at present using egcs1.1.1. we therefore
   separate the different implementations for the different dimensions using
   the preprocecssor and double check using an assertion in the function body.
*/
template <int dim>
template <typename number>
void
DoFCellAccessor<dim>::set_dof_values (const Vector<number> &local_values,
				      Vector<number>       &values) const {
  Assert (dim==1, ExcInternalError());
  
  Assert (dof_handler != 0, DoFAccessor<1>::ExcInvalidObject());
  Assert (&dof_handler->get_fe() != 0, DoFAccessor<1>::ExcInvalidObject());
  Assert (local_values.size() == dof_handler->get_fe().total_dofs,
	  DoFAccessor<1>::ExcVectorDoesNotMatch());
  Assert (values.size() == dof_handler->n_dofs(),
	  DoFAccessor<1>::ExcVectorDoesNotMatch());
  Assert (active(), ExcNotActive());
  
  const unsigned int dofs_per_vertex = dof_handler->get_fe().dofs_per_vertex,
		     dofs_per_line   = dof_handler->get_fe().dofs_per_line;
  typename vector<number>::const_iterator next_local_value=local_values.begin();
  for (unsigned int vertex=0; vertex<2; ++vertex)
    for (unsigned int d=0; d<dofs_per_vertex; ++d)
       values(vertex_dof_index(vertex,d)) = *next_local_value++;
  for (unsigned int d=0; d<dofs_per_line; ++d)
    values(dof_index(d)) = *next_local_value++;

  Assert (next_local_value == local_values.end(),
	  ::ExcInternalError());
};



#endif



#if deal_II_dimension == 2

template <>
TriaIterator<2, DoFObjectAccessor<1,2> >
DoFCellAccessor<2>::face (const unsigned int i) const
{
  return line(i);
};



/* this function looks like a template, but actually isn't one, it
   will only work for the dimension selected by the preprocessor
   guard above.
   
   the correct thing to do would be

   template <>
   template <typename number>
   void
   DoFCellAccessor<1>::get... (vector<number>...)

   but this does not work, at least not at present using egcs1.1.1. we therefore
   separate the different implementations for the different dimensions using
   the preprocecssor and double check using an assertion in the function body.
*/
template <int dim>
template <typename number>
void
DoFCellAccessor<dim>::get_dof_values (const Vector<number> &values,
				      Vector<number>       &local_values) const
{
  Assert (dim==2, ExcInternalError());
  
  Assert (dof_handler != 0, DoFAccessor<dim>::ExcInvalidObject());
  Assert (&dof_handler->get_fe() != 0, DoFAccessor<dim>::ExcInvalidObject());
  Assert (local_values.size() == dof_handler->get_fe().total_dofs,
	  DoFAccessor<dim>::ExcVectorDoesNotMatch());
  Assert (values.size() == dof_handler->n_dofs(),
	  DoFAccessor<dim>::ExcVectorDoesNotMatch());
  Assert (active(), ExcNotActive());
  
  const unsigned int dofs_per_vertex = dof_handler->get_fe().dofs_per_vertex,
		     dofs_per_line   = dof_handler->get_fe().dofs_per_line,
		     dofs_per_quad   = dof_handler->get_fe().dofs_per_quad;
  typename vector<number>::iterator next_local_value=local_values.begin();
  for (unsigned int vertex=0; vertex<4; ++vertex)
    for (unsigned int d=0; d<dofs_per_vertex; ++d)
      *next_local_value++ = values(vertex_dof_index(vertex,d));
  for (unsigned int line=0; line<4; ++line)
    for (unsigned int d=0; d<dofs_per_line; ++d)
      *next_local_value++ = values(this->line(line)->dof_index(d));
  for (unsigned int d=0; d<dofs_per_quad; ++d)
    *next_local_value++ = values(dof_index(d));

  Assert (next_local_value == local_values.end(),
	  ::ExcInternalError());
};



/* this function looks like a template, but actually isn't one, it
   will only work for the dimension selected by the preprocessor
   guard above.
   
   the correct thing to do would be

   template <>
   template <typename number>
   void
   DoFCellAccessor<1>::get... (vector<number>...)

   but this does not work, at least not at present using egcs1.1.1. we therefore
   separate the different implementations for the different dimensions using
   the preprocecssor and double check using an assertion in the function body.
*/
template <int dim>
template <typename number>
void
DoFCellAccessor<dim>::set_dof_values (const Vector<number> &local_values,
				      Vector<number>       &values) const {
  Assert (dim==2, ExcInternalError());
  
  Assert (dof_handler != 0, DoFAccessor<dim>::ExcInvalidObject());
  Assert (&dof_handler->get_fe() != 0, DoFAccessor<dim>::ExcInvalidObject());
  Assert (local_values.size() == dof_handler->get_fe().total_dofs,
	  DoFAccessor<dim>::ExcVectorDoesNotMatch());
  Assert (values.size() == dof_handler->n_dofs(),
	  DoFAccessor<dim>::ExcVectorDoesNotMatch());
  Assert (active(), ExcNotActive());
  
  const unsigned int dofs_per_vertex = dof_handler->get_fe().dofs_per_vertex,
		     dofs_per_line   = dof_handler->get_fe().dofs_per_line,
		     dofs_per_quad   = dof_handler->get_fe().dofs_per_quad;
  typename vector<number>::const_iterator next_local_value=local_values.begin();
  for (unsigned int vertex=0; vertex<4; ++vertex)
    for (unsigned int d=0; d<dofs_per_vertex; ++d)
      values(vertex_dof_index(vertex,d)) = *next_local_value++;
  for (unsigned int line=0; line<4; ++line)
    for (unsigned int d=0; d<dofs_per_line; ++d)
      values(this->line(line)->dof_index(d)) = *next_local_value++;
  for (unsigned int d=0; d<dofs_per_quad; ++d)
    values(dof_index(d)) = *next_local_value++;

  Assert (next_local_value == local_values.end(),
	  ::ExcInternalError());
};


#endif



#if deal_II_dimension == 3

template <>
TriaIterator<3, DoFObjectAccessor<2, 3> >
DoFCellAccessor<3>::face (const unsigned int i) const
{
  return quad(i);
};



/* this function looks like a template, but actually isn't one, it
   will only work for the dimension selected by the preprocessor
   guard above.
   
   the correct thing to do would be

   template <>
   template <typename number>
   void
   DoFCellAccessor<1>::get... (vector<number>...)

   but this does not work, at least not at present using egcs1.1.1. we therefore
   separate the different implementations for the different dimensions using
   the preprocecssor and double check using an assertion in the function body.
*/
template <int dim>
template <typename number>
void
DoFCellAccessor<dim>::get_dof_values (const Vector<number> &values,
				      Vector<number>       &local_values) const {
  Assert (dim==3, ExcInternalError());

  Assert (dof_handler != 0, DoFAccessor<3>::ExcInvalidObject());
  Assert (&dof_handler->get_fe() != 0, DoFAccessor<3>::ExcInvalidObject());
  Assert (local_values.size() == dof_handler->get_fe().total_dofs,
	  DoFAccessor<3>::ExcVectorDoesNotMatch());
  Assert (values.size() == dof_handler->n_dofs(),
	  DoFAccessor<3>::ExcVectorDoesNotMatch());
  Assert (active(), ExcNotActive());
  
  const unsigned int dofs_per_vertex = dof_handler->get_fe().dofs_per_vertex,
		     dofs_per_line   = dof_handler->get_fe().dofs_per_line,
		     dofs_per_quad   = dof_handler->get_fe().dofs_per_quad,
		     dofs_per_hex    = dof_handler->get_fe().dofs_per_hex;
  typename vector<number>::iterator next_local_value = local_values.begin();
  for (unsigned int vertex=0; vertex<GeometryInfo<dim>::vertices_per_cell; ++vertex)
    for (unsigned int d=0; d<dofs_per_vertex; ++d)
      *next_local_value++ = values(vertex_dof_index(vertex,d));
  for (unsigned int line=0; line<GeometryInfo<dim>::lines_per_cell; ++line)
    for (unsigned int d=0; d<dofs_per_line; ++d)
      *next_local_value++ = values(this->line(line)->dof_index(d));
  for (unsigned int quad=0; quad<GeometryInfo<dim>::quads_per_cell; ++quad)
    for (unsigned int d=0; d<dofs_per_quad; ++d)
      *next_local_value++ = values(this->quad(quad)->dof_index(d));
  for (unsigned int d=0; d<dofs_per_hex; ++d)
    *next_local_value++ = values(dof_index(d));

  Assert (next_local_value == local_values.end(),
	  ::ExcInternalError());
};



/* this function looks like a template, but actually isn't one, it
   will only work for the dimension selected by the preprocessor
   guard above.
   
   the correct thing to do would be

   template <>
   template <typename number>
   void
   DoFCellAccessor<1>::get... (vector<number>...)

   but this does not work, at least not at present using egcs1.1.1. we therefore
   separate the different implementations for the different dimensions using
   the preprocecssor and double check using an assertion in the function body.
*/
template <int dim>
template <typename number>
void
DoFCellAccessor<dim>::set_dof_values (const Vector<number> &local_values,
				      Vector<number>       &values) const {
  Assert (dim==3, ExcInternalError());

  Assert (dof_handler != 0, DoFAccessor<3>::ExcInvalidObject());
  Assert (&dof_handler->get_fe() != 0, DoFAccessor<3>::ExcInvalidObject());
  Assert (local_values.size() == dof_handler->get_fe().total_dofs,
	  DoFAccessor<3>::ExcVectorDoesNotMatch());
  Assert (values.size() == dof_handler->n_dofs(),
	  DoFAccessor<3>::ExcVectorDoesNotMatch());
  Assert (active(), ExcNotActive());
  
  const unsigned int dofs_per_vertex = dof_handler->get_fe().dofs_per_vertex,
		     dofs_per_line   = dof_handler->get_fe().dofs_per_line,
		     dofs_per_quad   = dof_handler->get_fe().dofs_per_quad,
		     dofs_per_hex    = dof_handler->get_fe().dofs_per_hex;
  typename vector<number>::const_iterator next_local_value=local_values.begin();
  for (unsigned int vertex=0; vertex<GeometryInfo<dim>::vertices_per_cell; ++vertex)
    for (unsigned int d=0; d<dofs_per_vertex; ++d)
      values(vertex_dof_index(vertex,d)) = *next_local_value++;
  for (unsigned int line=0; line<GeometryInfo<dim>::lines_per_cell; ++line)
    for (unsigned int d=0; d<dofs_per_line; ++d)
      values(this->line(line)->dof_index(d)) = *next_local_value++;
  for (unsigned int quad=0; quad<GeometryInfo<dim>::quads_per_cell; ++quad)
    for (unsigned int d=0; d<dofs_per_quad; ++d)
      values(this->quad(quad)->dof_index(d)) = *next_local_value++;
  for (unsigned int d=0; d<dofs_per_hex; ++d)
    values(dof_index(d)) = *next_local_value++;

  Assert (next_local_value == local_values.end(),
	  ::ExcInternalError());
};



#endif




template <int dim>
template <typename number>
void
DoFCellAccessor<dim>::get_interpolated_dof_values (const Vector<number> &values,
						   Vector<number>       &interpolated_values) const {
  const unsigned int total_dofs = dof_handler->get_fe().total_dofs;
  
  Assert (dof_handler != 0, DoFAccessor<dim>::ExcInvalidObject());
  Assert (&dof_handler->get_fe() != 0, DoFAccessor<dim>::ExcInvalidObject());
  Assert (interpolated_values.size() == total_dofs,
	  DoFAccessor<dim>::ExcVectorDoesNotMatch());
  Assert (values.size() == dof_handler->n_dofs(),
	  DoFAccessor<dim>::ExcVectorDoesNotMatch());

  if (!has_children())
				     // if this cell has no children: simply
				     // return the exact values on this cell
    get_dof_values (values, interpolated_values);
  else
				     // otherwise clobber them from the children
    {
      Vector<number> tmp1(total_dofs);
      Vector<number> tmp2(total_dofs);
      
      interpolated_values.clear ();

      const bool restriction_is_additive 
	= dof_handler->get_fe().restriction_is_additive;

      for (unsigned int child=0; child<GeometryInfo<dim>::children_per_cell;
	   ++child)
	{
					   // get the values from the present
					   // child, if necessary by
					   // interpolation itself
	  this->child(child)->get_interpolated_dof_values (values,
							   tmp1);
					   // interpolate these to the mother
					   // cell
	  dof_handler->get_fe().restrict(child).vmult (tmp2, tmp1);
	  
					   // now write those entries in tmp2
					   // which are != 0 into the output
					   // vector. Note that we may not
					   // add them up, since we would then
					   // end in adding up the contribution
					   // from nodes on boundaries of
					   // children more than once.
	  if (restriction_is_additive) 
	    {
	      for (unsigned int i=0; i<total_dofs; ++i)
		interpolated_values(i) += tmp2(i);
	    } 
	  else
	    {  
	      for (unsigned int i=0; i<total_dofs; ++i)
		if (tmp2(i) != 0)
		  interpolated_values(i) = tmp2(i);
	    };
	};
    };
};



template <int dim>
template <typename number>
void
DoFCellAccessor<dim>::set_dof_values_by_interpolation (const Vector<number> &local_values,
						       Vector<number>       &values) const {
  const unsigned int total_dofs = dof_handler->get_fe().total_dofs;
  
  Assert (dof_handler != 0, DoFAccessor<dim>::ExcInvalidObject());
  Assert (&dof_handler->get_fe() != 0, DoFAccessor<dim>::ExcInvalidObject());
  Assert (local_values.size() == total_dofs,
	  DoFAccessor<dim>::ExcVectorDoesNotMatch());
  Assert (values.size() == dof_handler->n_dofs(),
	  DoFAccessor<dim>::ExcVectorDoesNotMatch());

  if (!has_children())
				     // if this cell has no children: simply
				     // set the values on this cell
    set_dof_values (local_values, values);
  else
				     // otherwise distribute them to the children
    {
      Vector<number> tmp(total_dofs);
      
      for (unsigned int child=0; child<GeometryInfo<dim>::children_per_cell;
	   ++child)
	{
					   // prolong the given data
					   // to the present cell
	  dof_handler->get_fe().prolongate(child).vmult (tmp, local_values);
	  this->child(child)->set_dof_values_by_interpolation (tmp, values);
	};
    };
};







// explicit instantiations


// for double
template
void
DoFCellAccessor<deal_II_dimension>::get_dof_values (const Vector<double> &,
						    Vector<double>       &) const;
template
void
DoFCellAccessor<deal_II_dimension>::get_interpolated_dof_values (const Vector<double> &,
								 Vector<double>       &) const;
template
void
DoFCellAccessor<deal_II_dimension>::set_dof_values (const Vector<double> &,
						    Vector<double>       &) const;
template
void
DoFCellAccessor<deal_II_dimension>::set_dof_values_by_interpolation(const Vector<double> &,
								    Vector<double>       &) const;
// for float
template
void
DoFCellAccessor<deal_II_dimension>::get_dof_values (const Vector<float> &,
						    Vector<float>       &) const;
template
void
DoFCellAccessor<deal_II_dimension>::get_interpolated_dof_values (const Vector<float> &,
								 Vector<float>       &) const;
template
void
DoFCellAccessor<deal_II_dimension>::set_dof_values (const Vector<float> &,
						    Vector<float>       &) const;
template
void
DoFCellAccessor<deal_II_dimension>::set_dof_values_by_interpolation(const Vector<float> &,
								    Vector<float>       &) const;



#if deal_II_dimension == 1
template class DoFObjectAccessor<1, 1>;
template class DoFCellAccessor<1>;

template class TriaRawIterator<1,DoFCellAccessor<1> >;
template class TriaIterator<1,DoFCellAccessor<1> >;
template class TriaActiveIterator<1,DoFCellAccessor<1> >;
#endif

#if deal_II_dimension == 2
template class DoFObjectAccessor<1, 2>;
template class DoFObjectAccessor<2, 2>;
template class DoFCellAccessor<2>;

template class TriaRawIterator   <2,DoFObjectAccessor<1, 2> >;
template class TriaIterator      <2,DoFObjectAccessor<1, 2> >;
template class TriaActiveIterator<2,DoFObjectAccessor<1, 2> >;
template class TriaRawIterator   <2,DoFCellAccessor<2> >;
template class TriaIterator      <2,DoFCellAccessor<2> >;
template class TriaActiveIterator<2,DoFCellAccessor<2> >;
#endif



#if deal_II_dimension == 3
template class DoFObjectAccessor<1, 3>;
template class DoFObjectAccessor<2, 3>;
template class DoFCellAccessor<3>;

template class TriaRawIterator   <3,DoFObjectAccessor<1, 3> >;
template class TriaIterator      <3,DoFObjectAccessor<1, 3> >;
template class TriaActiveIterator<3,DoFObjectAccessor<1, 3> >;
template class TriaRawIterator   <3,DoFObjectAccessor<2, 3> >;
template class TriaIterator      <3,DoFObjectAccessor<2, 3> >;
template class TriaActiveIterator<3,DoFObjectAccessor<2, 3> >;
template class TriaRawIterator   <3,DoFCellAccessor<3> >;
template class TriaIterator      <3,DoFCellAccessor<3> >;
template class TriaActiveIterator<3,DoFCellAccessor<3> >;
#endif


