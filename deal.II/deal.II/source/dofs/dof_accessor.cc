//----------------------------  dof_accessor.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  dof_accessor.cc  ---------------------------


#include <lac/vector.h>
#include <lac/block_vector.h>
#include <lac/sparse_matrix.h>

#include <dofs/dof_accessor.h>
#include <dofs/dof_accessor.templates.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_levels.h>
#include <grid/tria_iterator.h>
#include <grid/tria_iterator.templates.h>
#include <fe/fe.h>

#include <vector>


/*------------------------- Functions: DoFLineAccessor -----------------------*/


template <int dim>
void DoFObjectAccessor<1, dim>::set_dof_index (const unsigned int i,
					       const unsigned int index) const
{
  Assert (this->dof_handler != 0, DoFAccessor<dim>::ExcInvalidObject());
				   // make sure a FE has been selected
				   // and enough room was reserved
  Assert (this->dof_handler->selected_fe != 0, DoFAccessor<dim>::ExcInvalidObject());
  Assert (i<this->dof_handler->selected_fe->dofs_per_line,
	  ExcIndexRange (i, 0, this->dof_handler->selected_fe->dofs_per_line));

  this->dof_handler->levels[this->present_level]
    ->line_dofs[this->present_index*this->dof_handler->selected_fe->dofs_per_line+i] = index;
};



template <int dim>
void DoFObjectAccessor<1, dim>::set_vertex_dof_index (const unsigned int vertex,
						      const unsigned int i,
						      const unsigned int index) const

{
				   // since the exception classes are
				   // from a template dependent base
				   // class, we have to fully qualify
				   // them. to work around more
				   // trouble, typedef the template
				   // dependent base class to a
				   // non-template dependent name and
				   // use that to specify the
				   // qualified exception names
  typedef DoFAccessor<dim> BaseClass;
  
  Assert (this->dof_handler != 0,
	  typename BaseClass::ExcInvalidObject());
  Assert (this->dof_handler->selected_fe != 0,
	  typename BaseClass::ExcInvalidObject());
  Assert (vertex<2, ExcIndexRange (i,0,2));
  Assert (i<this->dof_handler->selected_fe->dofs_per_vertex,
	  ExcIndexRange (i, 0, this->dof_handler->selected_fe->dofs_per_vertex));

  const unsigned int dof_number = (this->vertex_index(vertex) *
				   this->dof_handler->selected_fe->dofs_per_vertex +
				   i);
  this->dof_handler->vertex_dofs[dof_number] = index;
};



template <int dim>
void DoFObjectAccessor<1, dim>::
distribute_local_to_global (const Vector<double> &local_source,
			    Vector<double>       &global_destination) const
{
				   // since the exception classes are
				   // from a template dependent base
				   // class, we have to fully qualify
				   // them. to work around more
				   // trouble, typedef the template
				   // dependent base class to a
				   // non-template dependent name and
				   // use that to specify the
				   // qualified exception names
  typedef DoFAccessor<dim> BaseClass;
  
  Assert (this->dof_handler != 0,
	  typename BaseClass::ExcInvalidObject());
  Assert (this->dof_handler->selected_fe != 0,
	  typename BaseClass::ExcInvalidObject());
  Assert (local_source.size() == (2*this->dof_handler->get_fe().dofs_per_vertex +
				  this->dof_handler->get_fe().dofs_per_line),
	  typename BaseClass::ExcVectorDoesNotMatch());
  Assert (this->dof_handler->n_dofs() == global_destination.size(),
	  typename BaseClass::ExcVectorDoesNotMatch());

  const unsigned int n_dofs = local_source.size();
  
				   // get indices of dofs
  std::vector<unsigned int> dofs (n_dofs);
  get_dof_indices (dofs);
  
				   // distribute cell vector
  for (unsigned int j=0; j<n_dofs; ++j)
    global_destination(dofs[j]) += local_source(j);
};



template <int dim>
void DoFObjectAccessor<1, dim>::
distribute_local_to_global (const FullMatrix<double> &local_source,
			    SparseMatrix<double>     &global_destination) const
{
				   // since the exception classes are
				   // from a template dependent base
				   // class, we have to fully qualify
				   // them. to work around more
				   // trouble, typedef the template
				   // dependent base class to a
				   // non-template dependent name and
				   // use that to specify the
				   // qualified exception names
  typedef DoFAccessor<dim> BaseClass;
  
  Assert (this->dof_handler != 0,
	  typename BaseClass::ExcInvalidObject());
  Assert (this->dof_handler->selected_fe != 0,
	  typename BaseClass::ExcInvalidObject());
  Assert (local_source.m() == (2*this->dof_handler->get_fe().dofs_per_vertex +
			       this->dof_handler->get_fe().dofs_per_line),
	  typename BaseClass::ExcVectorDoesNotMatch());
  Assert (local_source.m() == local_source.n(),
	  typename BaseClass::ExcMatrixDoesNotMatch());
  Assert (this->dof_handler->n_dofs() == global_destination.m(),
	  typename BaseClass::ExcMatrixDoesNotMatch());
  Assert (global_destination.m() == global_destination.n(),
	  typename BaseClass::ExcMatrixDoesNotMatch());

  const unsigned int n_dofs = local_source.m();

				   // get indices of dofs
  std::vector<unsigned int> dofs (n_dofs);
  get_dof_indices (dofs);
  
				   // distribute cell matrix
  for (unsigned int i=0; i<n_dofs; ++i)
    for (unsigned int j=0; j<n_dofs; ++j)
      global_destination.add(dofs[i], dofs[j], local_source(i,j));
};



template <int dim>
template <class InputVector, typename number>
void
DoFObjectAccessor<1,dim>::get_dof_values (const InputVector &values,
					  Vector<number>    &local_values) const
{
  Assert (dim==1, ExcInternalError());
  
  Assert (this->dof_handler != 0, DoFAccessor<1>::ExcInvalidObject());
  Assert (&this->dof_handler->get_fe() != 0, DoFAccessor<1>::ExcInvalidObject());
  Assert (local_values.size() == this->dof_handler->get_fe().dofs_per_cell,
	  DoFAccessor<1>::ExcVectorDoesNotMatch());
  Assert (values.size() == this->dof_handler->n_dofs(),
	  DoFAccessor<1>::ExcVectorDoesNotMatch());
  Assert (this->has_children() == false, DoFAccessor<1>::ExcNotActive());
  
  const unsigned int dofs_per_vertex = this->dof_handler->get_fe().dofs_per_vertex,
		     dofs_per_line   = this->dof_handler->get_fe().dofs_per_line;
  typename Vector<number>::iterator next_local_value=local_values.begin();
  for (unsigned int vertex=0; vertex<2; ++vertex)
    for (unsigned int d=0; d<dofs_per_vertex; ++d)
      *next_local_value++ = values(vertex_dof_index(vertex,d));
  for (unsigned int d=0; d<dofs_per_line; ++d)
    *next_local_value++ = values(dof_index(d));

  Assert (next_local_value == local_values.end(),
	  ExcInternalError());
};



template <int dim>
template <class OutputVector, typename number>
void
DoFObjectAccessor<1,dim>::set_dof_values (const Vector<number> &local_values,
					  OutputVector         &values) const
{
  Assert (dim==1, ExcInternalError());
  
  Assert (this->dof_handler != 0, DoFAccessor<1>::ExcInvalidObject());
  Assert (&this->dof_handler->get_fe() != 0, DoFAccessor<1>::ExcInvalidObject());
  Assert (local_values.size() == this->dof_handler->get_fe().dofs_per_cell,
	  DoFAccessor<1>::ExcVectorDoesNotMatch());
  Assert (values.size() == this->dof_handler->n_dofs(),
	  DoFAccessor<1>::ExcVectorDoesNotMatch());
  Assert (this->has_children() == false, DoFAccessor<1>::ExcNotActive());
  
  const unsigned int dofs_per_vertex = this->dof_handler->get_fe().dofs_per_vertex,
		     dofs_per_line   = this->dof_handler->get_fe().dofs_per_line;
  typename Vector<number>::const_iterator next_local_value=local_values.begin();
  for (unsigned int vertex=0; vertex<2; ++vertex)
    for (unsigned int d=0; d<dofs_per_vertex; ++d)
      values(vertex_dof_index(vertex,d)) = *next_local_value++;
  for (unsigned int d=0; d<dofs_per_line; ++d)
    values(dof_index(d)) = *next_local_value++;

  Assert (next_local_value == local_values.end(),
	  ExcInternalError());
};


/*------------------------- Functions: DoFQuadAccessor -----------------------*/


template <int dim>
void DoFObjectAccessor<2, dim>::set_dof_index (const unsigned int i,
					       const unsigned int index) const
{
  Assert (this->dof_handler != 0,
	  typename DoFAccessor<dim>::ExcInvalidObject());
				   // make sure a FE has been selected
				   // and enough room was reserved
  Assert (this->dof_handler->selected_fe != 0,
	  typename DoFAccessor<dim>::ExcInvalidObject());
  Assert (i<this->dof_handler->selected_fe->dofs_per_quad,
	  ExcIndexRange (i, 0, this->dof_handler->selected_fe->dofs_per_quad));

  this->dof_handler->levels[this->present_level]
    ->quad_dofs[this->present_index*this->dof_handler->selected_fe->dofs_per_quad+i] = index;
};



template <int dim>
void
DoFObjectAccessor<2, dim>::set_vertex_dof_index (const unsigned int vertex,
						 const unsigned int i,
						 const unsigned int index) const
{
  Assert (this->dof_handler != 0,
	  typename DoFAccessor<dim>::ExcInvalidObject());
  Assert (this->dof_handler->selected_fe != 0,
	  typename DoFAccessor<dim>::ExcInvalidObject());
  Assert (vertex<4, ExcIndexRange (i,0,4));
  Assert (i<this->dof_handler->selected_fe->dofs_per_vertex,
	  ExcIndexRange (i, 0, this->dof_handler->selected_fe->dofs_per_vertex));

  const unsigned int dof_number = (this->vertex_index(vertex) *
				   this->dof_handler->selected_fe->dofs_per_vertex +
				   i);
  this->dof_handler->vertex_dofs[dof_number] = index;
};



template <int dim>
void DoFObjectAccessor<2, dim>::
distribute_local_to_global (const Vector<double> &local_source,
			    Vector<double>       &global_destination) const {
  Assert (this->dof_handler != 0,
	  typename DoFAccessor<dim>::ExcInvalidObject());
  Assert (this->dof_handler->selected_fe != 0,
	  typename DoFAccessor<dim>::ExcInvalidObject());
  Assert (local_source.size() == (4*this->dof_handler->get_fe().dofs_per_vertex +
				  4*this->dof_handler->get_fe().dofs_per_line +
				  this->dof_handler->get_fe().dofs_per_quad),
	  typename DoFAccessor<dim>::ExcVectorDoesNotMatch());
  Assert (this->dof_handler->n_dofs() == global_destination.size(),
	  typename DoFAccessor<dim>::ExcVectorDoesNotMatch());

  const unsigned int n_dofs = local_source.size();

				   // get indices of dofs
  std::vector<unsigned int> dofs (n_dofs);
  get_dof_indices (dofs);
  
				   // distribute cell vector
  for (unsigned int j=0; j<n_dofs; ++j)
    global_destination(dofs[j]) += local_source(j);
};



template <int dim>
void DoFObjectAccessor<2, dim>::
distribute_local_to_global (const FullMatrix<double> &local_source,
			    SparseMatrix<double>     &global_destination) const {
  Assert (this->dof_handler != 0,
	  typename DoFAccessor<dim>::ExcInvalidObject());
  Assert (this->dof_handler->selected_fe != 0,
	  typename DoFAccessor<dim>::ExcInvalidObject());
  Assert (local_source.m() == (4*this->dof_handler->get_fe().dofs_per_vertex +
			       4*this->dof_handler->get_fe().dofs_per_line +
			       this->dof_handler->get_fe().dofs_per_quad),
	  typename DoFAccessor<dim>::ExcMatrixDoesNotMatch());
  Assert (local_source.m() == local_source.n(),
	  typename DoFAccessor<dim>::ExcMatrixDoesNotMatch());
  Assert (this->dof_handler->n_dofs() == global_destination.m(),
	  typename DoFAccessor<dim>::ExcMatrixDoesNotMatch());
  Assert (global_destination.m() == global_destination.n(),
	  typename DoFAccessor<dim>::ExcMatrixDoesNotMatch());
  
  const unsigned int n_dofs = local_source.m();

				   // get indices of dofs
  std::vector<unsigned int> dofs (n_dofs);
  get_dof_indices (dofs);
  
				   // distribute cell matrix
  for (unsigned int i=0; i<n_dofs; ++i)
    for (unsigned int j=0; j<n_dofs; ++j)
      global_destination.add(dofs[i], dofs[j], local_source(i,j));
};



template <int dim>
template <class InputVector, typename number>
void
DoFObjectAccessor<2,dim>::get_dof_values (const InputVector &values,
					  Vector<number>    &local_values) const
{
  Assert (dim==2, ExcInternalError());
  
  Assert (this->dof_handler != 0,
	  typename DoFAccessor<dim>::ExcInvalidObject());
  Assert (&this->dof_handler->get_fe() != 0,
	  typename DoFAccessor<dim>::ExcInvalidObject());
  Assert (local_values.size() == this->dof_handler->get_fe().dofs_per_cell,
	  typename DoFAccessor<dim>::ExcVectorDoesNotMatch());
  Assert (values.size() == this->dof_handler->n_dofs(),
	  typename DoFAccessor<dim>::ExcVectorDoesNotMatch());
  Assert (this->has_children() == false, typename DoFAccessor<dim>::
	  ExcNotActive());
  
  const unsigned int dofs_per_vertex = this->dof_handler->get_fe().dofs_per_vertex,
		     dofs_per_line   = this->dof_handler->get_fe().dofs_per_line,
		     dofs_per_quad   = this->dof_handler->get_fe().dofs_per_quad;
  typename Vector<number>::iterator next_local_value=local_values.begin();
  for (unsigned int vertex=0; vertex<4; ++vertex)
    for (unsigned int d=0; d<dofs_per_vertex; ++d)
      *next_local_value++ = values(vertex_dof_index(vertex,d));
  for (unsigned int line=0; line<4; ++line)
    for (unsigned int d=0; d<dofs_per_line; ++d)
      *next_local_value++ = values(this->line(line)->dof_index(d));
  for (unsigned int d=0; d<dofs_per_quad; ++d)
    *next_local_value++ = values(dof_index(d));

  Assert (next_local_value == local_values.end(),
	  ExcInternalError());
};



template <int dim>
template <class OutputVector, typename number>
void
DoFObjectAccessor<2,dim>::set_dof_values (const Vector<number> &local_values,
					  OutputVector         &values) const
{
  Assert (dim==2, ExcInternalError());
  
  Assert (this->dof_handler != 0,
	  typename DoFAccessor<dim>::ExcInvalidObject());
  Assert (&this->dof_handler->get_fe() != 0,
	  typename DoFAccessor<dim>::ExcInvalidObject());
  Assert (local_values.size() == this->dof_handler->get_fe().dofs_per_cell,
	  typename DoFAccessor<dim>::ExcVectorDoesNotMatch());
  Assert (values.size() == this->dof_handler->n_dofs(),
	  typename DoFAccessor<dim>::ExcVectorDoesNotMatch());
  Assert (this->has_children() == false,
	  typename DoFAccessor<dim>::ExcNotActive());
  
  const unsigned int dofs_per_vertex = this->dof_handler->get_fe().dofs_per_vertex,
		     dofs_per_line   = this->dof_handler->get_fe().dofs_per_line,
		     dofs_per_quad   = this->dof_handler->get_fe().dofs_per_quad;
  typename Vector<number>::const_iterator next_local_value=local_values.begin();
  for (unsigned int vertex=0; vertex<4; ++vertex)
    for (unsigned int d=0; d<dofs_per_vertex; ++d)
      values(vertex_dof_index(vertex,d)) = *next_local_value++;
  for (unsigned int line=0; line<4; ++line)
    for (unsigned int d=0; d<dofs_per_line; ++d)
      values(this->line(line)->dof_index(d)) = *next_local_value++;
  for (unsigned int d=0; d<dofs_per_quad; ++d)
    values(dof_index(d)) = *next_local_value++;

  Assert (next_local_value == local_values.end(),
	  ExcInternalError());
};


/*------------------------- Functions: DoFHexAccessor -----------------------*/


template <int dim>
void DoFObjectAccessor<3, dim>::set_dof_index (const unsigned int i,
					       const unsigned int index) const
{
  Assert (this->dof_handler != 0,
	  typename DoFAccessor<dim>::ExcInvalidObject());
				   // make sure a FE has been selected
				   // and enough room was reserved
  Assert (this->dof_handler->selected_fe != 0,
	  typename DoFAccessor<dim>::ExcInvalidObject());
  Assert (i<this->dof_handler->selected_fe->dofs_per_hex,
	  ExcIndexRange (i, 0, this->dof_handler->selected_fe->dofs_per_hex));

  this->dof_handler->levels[this->present_level]
    ->hex_dofs[this->present_index*this->dof_handler->selected_fe->dofs_per_hex+i] = index;
};



template <int dim>
void DoFObjectAccessor<3, dim>::set_vertex_dof_index (const unsigned int vertex,
						      const unsigned int i,
						      const unsigned int index) const
{
  Assert (this->dof_handler != 0,
	  typename DoFAccessor<dim>::ExcInvalidObject());
  Assert (this->dof_handler->selected_fe != 0,
	  typename DoFAccessor<dim>::ExcInvalidObject());
  Assert (vertex<8,
	  ExcIndexRange (i,0,8));
  Assert (i<this->dof_handler->selected_fe->dofs_per_vertex,
	  ExcIndexRange (i, 0, this->dof_handler->selected_fe->dofs_per_vertex));

  const unsigned int dof_number = (this->vertex_index(vertex) *
				   this->dof_handler->selected_fe->dofs_per_vertex +
				   i);
  this->dof_handler->vertex_dofs[dof_number] = index;
};



template <int dim>
void DoFObjectAccessor<3, dim>::
distribute_local_to_global (const Vector<double> &local_source,
			    Vector<double>       &global_destination) const
{
  Assert (this->dof_handler != 0,
	  typename DoFAccessor<dim>::ExcInvalidObject());
  Assert (this->dof_handler->selected_fe != 0,
	  typename DoFAccessor<dim>::ExcInvalidObject());
  Assert (local_source.size() == (8*this->dof_handler->get_fe().dofs_per_vertex +
				  12*this->dof_handler->get_fe().dofs_per_line +
				  6*this->dof_handler->get_fe().dofs_per_quad +
				  this->dof_handler->get_fe().dofs_per_hex),
	  typename DoFAccessor<dim>::ExcVectorDoesNotMatch());
  Assert (this->dof_handler->n_dofs() == global_destination.size(),
	  typename DoFAccessor<dim>::ExcVectorDoesNotMatch());

  const unsigned int n_dofs = local_source.size();

				   // get indices of dofs
  std::vector<unsigned int> dofs (n_dofs);
  get_dof_indices (dofs);
  
				   // distribute cell vector
  for (unsigned int j=0; j<n_dofs; ++j)
    global_destination(dofs[j]) += local_source(j);
};



template <int dim>
void DoFObjectAccessor<3, dim>::
distribute_local_to_global (const FullMatrix<double> &local_source,
			    SparseMatrix<double>     &global_destination) const
{
  Assert (this->dof_handler != 0,
	  typename DoFAccessor<dim>::ExcInvalidObject());
  Assert (this->dof_handler->selected_fe != 0,
	  typename DoFAccessor<dim>::ExcInvalidObject());
  Assert (local_source.m() == (8*this->dof_handler->get_fe().dofs_per_vertex +
			       12*this->dof_handler->get_fe().dofs_per_line +
			       6*this->dof_handler->get_fe().dofs_per_quad +
			       this->dof_handler->get_fe().dofs_per_hex),
	  typename DoFAccessor<dim>::ExcMatrixDoesNotMatch());
  Assert (local_source.m() == local_source.n(),
	  typename DoFAccessor<dim>::ExcMatrixDoesNotMatch());
  Assert (this->dof_handler->n_dofs() == global_destination.m(),
	  typename DoFAccessor<dim>::ExcMatrixDoesNotMatch());
  Assert (global_destination.m() == global_destination.n(),
	  typename DoFAccessor<dim>::ExcMatrixDoesNotMatch());
  
  const unsigned int n_dofs = local_source.m();

				   // get indices of dofs
  std::vector<unsigned int> dofs (n_dofs);
  get_dof_indices (dofs);
  
				   // distribute cell matrix
  for (unsigned int i=0; i<n_dofs; ++i)
    for (unsigned int j=0; j<n_dofs; ++j)
      global_destination.add(dofs[i], dofs[j], local_source(i,j));
};



template <int dim>
template <class InputVector, typename number>
void
DoFObjectAccessor<3,dim>::get_dof_values (const InputVector &values,
					  Vector<number>    &local_values) const
{
  Assert (dim==3, ExcInternalError());

  Assert (this->dof_handler != 0, DoFAccessor<3>::ExcInvalidObject());
  Assert (&this->dof_handler->get_fe() != 0, DoFAccessor<3>::ExcInvalidObject());
  Assert (local_values.size() == this->dof_handler->get_fe().dofs_per_cell,
	  DoFAccessor<3>::ExcVectorDoesNotMatch());
  Assert (values.size() == this->dof_handler->n_dofs(),
	  DoFAccessor<3>::ExcVectorDoesNotMatch());
  Assert (this->has_children() == false,
	  DoFAccessor<3>::ExcNotActive());
  
  const unsigned int dofs_per_vertex = this->dof_handler->get_fe().dofs_per_vertex,
		     dofs_per_line   = this->dof_handler->get_fe().dofs_per_line,
		     dofs_per_quad   = this->dof_handler->get_fe().dofs_per_quad,
		     dofs_per_hex    = this->dof_handler->get_fe().dofs_per_hex;
  typename Vector<number>::iterator next_local_value = local_values.begin();
  for (unsigned int vertex=0; vertex<GeometryInfo<3>::vertices_per_cell; ++vertex)
    for (unsigned int d=0; d<dofs_per_vertex; ++d)
      *next_local_value++ = values(vertex_dof_index(vertex,d));
  for (unsigned int line=0; line<GeometryInfo<3>::lines_per_cell; ++line)
    for (unsigned int d=0; d<dofs_per_line; ++d)
      *next_local_value++ = values(this->line(line)->dof_index(d));
  for (unsigned int quad=0; quad<GeometryInfo<3>::quads_per_cell; ++quad)
    for (unsigned int d=0; d<dofs_per_quad; ++d)
      *next_local_value++ = values(this->quad(quad)->dof_index(d));
  for (unsigned int d=0; d<dofs_per_hex; ++d)
    *next_local_value++ = values(dof_index(d));

  Assert (next_local_value == local_values.end(),
	  ExcInternalError());
};



template <int dim>
template <class OutputVector, typename number>
void
DoFObjectAccessor<3,dim>::set_dof_values (const Vector<number> &local_values,
					  OutputVector         &values) const
{
  Assert (dim==3, ExcInternalError());

  Assert (this->dof_handler != 0,
	  DoFAccessor<3>::ExcInvalidObject());
  Assert (&this->dof_handler->get_fe() != 0,
	  DoFAccessor<3>::ExcInvalidObject());
  Assert (local_values.size() == this->dof_handler->get_fe().dofs_per_cell,
	  DoFAccessor<3>::ExcVectorDoesNotMatch());
  Assert (values.size() == this->dof_handler->n_dofs(),
	  DoFAccessor<3>::ExcVectorDoesNotMatch());
  Assert (this->has_children() == false,
	  DoFAccessor<3>::ExcNotActive());
  
  const unsigned int dofs_per_vertex = this->dof_handler->get_fe().dofs_per_vertex,
		     dofs_per_line   = this->dof_handler->get_fe().dofs_per_line,
		     dofs_per_quad   = this->dof_handler->get_fe().dofs_per_quad,
		     dofs_per_hex    = this->dof_handler->get_fe().dofs_per_hex;
  typename Vector<number>::const_iterator next_local_value=local_values.begin();
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
	  ExcInternalError());
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

#endif


#if deal_II_dimension == 2

template <>
TriaIterator<2, DoFObjectAccessor<1,2> >
DoFCellAccessor<2>::face (const unsigned int i) const
{
  return this->line(i);
};

#endif


#if deal_II_dimension == 3

template <>
TriaIterator<3, DoFObjectAccessor<2, 3> >
DoFCellAccessor<3>::face (const unsigned int i) const
{
  return this->quad(i);
};

#endif


template <int dim>
template <class InputVector, typename number>
void
DoFCellAccessor<dim>::get_interpolated_dof_values (const InputVector &values,
						   Vector<number>    &interpolated_values) const
{
  const FiniteElement<dim> &fe            = this->dof_handler->get_fe();
  const unsigned int        dofs_per_cell = fe.dofs_per_cell;
  
  
  Assert (this->dof_handler != 0, DoFAccessor<dim>::ExcInvalidObject());
  Assert (&fe != 0, DoFAccessor<dim>::ExcInvalidObject());
  Assert (interpolated_values.size() == dofs_per_cell,
	  DoFAccessor<dim>::ExcVectorDoesNotMatch());
  Assert (values.size() == this->dof_handler->n_dofs(),
	  DoFAccessor<dim>::ExcVectorDoesNotMatch());

  if (!this->has_children())
				     // if this cell has no children: simply
				     // return the exact values on this cell
    get_dof_values (values, interpolated_values);
  else
				     // otherwise clobber them from the children
    {
      Vector<number> tmp1(dofs_per_cell);
      Vector<number> tmp2(dofs_per_cell);
      
      interpolated_values.clear ();

                                       // later on we will have to
                                       // push the values interpolated
                                       // from the child to the mother
                                       // cell into the output
                                       // vector. unfortunately, there
                                       // are two types of elements:
                                       // ones where you add up the
                                       // contributions from the
                                       // different child cells, and
                                       // ones where you overwrite.
                                       //
                                       // an example for the first is
                                       // piecewise constant (and
                                       // discontinuous) elements,
                                       // where we build the value on
                                       // the coarse cell by averaging
                                       // the values from the cell
                                       // (i.e. by adding up a
                                       // fraction of the values of
                                       // their values)
                                       //
                                       // and example for the latter
                                       // are the usual continuous
                                       // elements. the value on a
                                       // vertex of a coarse cell must
                                       // there be the same,
                                       // irrespective of the adjacent
                                       // cell we are presently on. so
                                       // we always overwrite. in
                                       // fact, we must, since we
                                       // cannot know in advance how
                                       // many neighbors there will
                                       // be, so there is no way to
                                       // compute the average with
                                       // fixed factors
                                       //
                                       // so we have to find out to
                                       // which type this element
                                       // belongs. the difficulty is:
                                       // the finite element may be a
                                       // composed one, so we can only
                                       // hope to do this for each
                                       // shape function
                                       // individually. to avoid doing
                                       // this over and over again, we
                                       // do this once now and cache
                                       // the results
      std::vector<bool> restriction_is_additive (dofs_per_cell);
      for (unsigned int i=0; i<dofs_per_cell; ++i)
        restriction_is_additive[i] = fe.restriction_is_additive(i);
      
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
	  fe.restrict(child).vmult (tmp2, tmp1);

                                           // and add up or set them
                                           // in the output vector
	  for (unsigned int i=0; i<dofs_per_cell; ++i)
            if (restriction_is_additive[i]) 
              interpolated_values(i) += tmp2(i);
            else
              if (tmp2(i) != 0)
                interpolated_values(i) = tmp2(i);
	};
    };
};



template <int dim>
template <class OutputVector, typename number>
void
DoFCellAccessor<dim>::set_dof_values_by_interpolation (const Vector<number> &local_values,
						       OutputVector         &values) const
{
  const unsigned int dofs_per_cell = this->dof_handler->get_fe().dofs_per_cell;
  
  Assert (this->dof_handler != 0, DoFAccessor<dim>::ExcInvalidObject());
  Assert (&this->dof_handler->get_fe() != 0, DoFAccessor<dim>::ExcInvalidObject());
  Assert (local_values.size() == dofs_per_cell,
	  DoFAccessor<dim>::ExcVectorDoesNotMatch());
  Assert (values.size() == this->dof_handler->n_dofs(),
	  DoFAccessor<dim>::ExcVectorDoesNotMatch());

  if (!this->has_children())
                                     // if this cell has no children: simply
				     // set the values on this cell
    set_dof_values (local_values, values);
  else
				     // otherwise distribute them to the children
    {
      Vector<number> tmp(dofs_per_cell);
      
      for (unsigned int child=0; child<GeometryInfo<dim>::children_per_cell;
	   ++child)
	{
					   // prolong the given data
					   // to the present cell
	  this->dof_handler->get_fe().prolongate(child).vmult (tmp, local_values);
	  this->child(child)->set_dof_values_by_interpolation (tmp, values);
	};
    };
};



// --------------------------------------------------------------------------
// explicit instantiations


template
void
DoFObjectAccessor<1,deal_II_dimension>::get_dof_values<Vector<double>,double>
(const Vector<double>&, Vector<double>&) const;
template
void
DoFObjectAccessor<1,deal_II_dimension>::set_dof_values<Vector<double>,double>
(const Vector<double>&, Vector<double>&) const;
template
void
DoFObjectAccessor<1,deal_II_dimension>::get_dof_values<Vector<double>,float>
(const Vector<double>&, Vector<float>&) const;
template
void
DoFObjectAccessor<1,deal_II_dimension>::set_dof_values<Vector<float>,double>
(const Vector<double>&, Vector<float>&) const;
template
void
DoFObjectAccessor<1,deal_II_dimension>::get_dof_values<Vector<float>,double>
(const Vector<float>&, Vector<double>&) const;
template
void
DoFObjectAccessor<1,deal_II_dimension>::set_dof_values<Vector<double>,float>
(const Vector<float>&, Vector<double>&) const;
template
void
DoFObjectAccessor<1,deal_II_dimension>::get_dof_values<Vector<float>,float>
(const Vector<float>&, Vector<float>&) const;
template
void
DoFObjectAccessor<1,deal_II_dimension>::set_dof_values<Vector<float>,float>
(const Vector<float>&, Vector<float>&) const;


// for block vector
template
void
DoFObjectAccessor<1,deal_II_dimension>::get_dof_values<BlockVector<double>,double>
(const BlockVector<double> &, Vector<double>&) const;
template
void
DoFObjectAccessor<1,deal_II_dimension>::set_dof_values<BlockVector<double>,double>
(const Vector<double>&, BlockVector<double>&) const;
template
void
DoFObjectAccessor<1,deal_II_dimension>::get_dof_values<BlockVector<double>,float>
(const BlockVector<double> &, Vector<float>&) const;
template
void
DoFObjectAccessor<1,deal_II_dimension>::set_dof_values<BlockVector<float>,double>
(const Vector<double>&, BlockVector<float>&) const;
template
void
DoFObjectAccessor<1,deal_II_dimension>::get_dof_values<BlockVector<float>,double>
(const BlockVector<float> &, Vector<double>&) const;
template
void
DoFObjectAccessor<1,deal_II_dimension>::set_dof_values<BlockVector<double>,float>
(const Vector<float>&, BlockVector<double>&) const;
template
void
DoFObjectAccessor<1,deal_II_dimension>::get_dof_values<BlockVector<float>,float>
(const BlockVector<float>&, Vector<float>&) const;
template
void
DoFObjectAccessor<1,deal_II_dimension>::set_dof_values<BlockVector<float>,float>
(const Vector<float>&, BlockVector<float>&) const;

#if deal_II_dimension >= 2
template
void
DoFObjectAccessor<2,deal_II_dimension>::get_dof_values<Vector<double>,double>
(const Vector<double>&, Vector<double>&) const;
template
void
DoFObjectAccessor<2,deal_II_dimension>::set_dof_values<Vector<double>,double>
(const Vector<double>&, Vector<double>&) const;
template
void
DoFObjectAccessor<2,deal_II_dimension>::get_dof_values<Vector<double>,float>
(const Vector<double>&, Vector<float>&) const;
template
void
DoFObjectAccessor<2,deal_II_dimension>::set_dof_values<Vector<float>,double>
(const Vector<double>&, Vector<float>&) const;
template
void
DoFObjectAccessor<2,deal_II_dimension>::get_dof_values<Vector<float>,double>
(const Vector<float>&, Vector<double>&) const;
template
void
DoFObjectAccessor<2,deal_II_dimension>::set_dof_values<Vector<double>,float>
(const Vector<float>&, Vector<double>&) const;
template
void
DoFObjectAccessor<2,deal_II_dimension>::get_dof_values<Vector<float>,float>
(const Vector<float>&, Vector<float>&) const;
template
void
DoFObjectAccessor<2,deal_II_dimension>::set_dof_values<Vector<float>,float>
(const Vector<float>&, Vector<float>&) const;


// for block vector
template
void
DoFObjectAccessor<2,deal_II_dimension>::get_dof_values<BlockVector<double>,double>
(const BlockVector<double>&, Vector<double>&) const;
template
void
DoFObjectAccessor<2,deal_II_dimension>::set_dof_values<BlockVector<double>,double>
(const Vector<double>&, BlockVector<double>&) const;
template
void
DoFObjectAccessor<2,deal_II_dimension>::get_dof_values<BlockVector<double>,float>
(const BlockVector<double>&, Vector<float>&) const;
template
void
DoFObjectAccessor<2,deal_II_dimension>::set_dof_values<BlockVector<float>,double>
(const Vector<double>&, BlockVector<float>&) const;
template
void
DoFObjectAccessor<2,deal_II_dimension>::get_dof_values<BlockVector<float>,double>
(const BlockVector<float>&, Vector<double>&) const;
template
void
DoFObjectAccessor<2,deal_II_dimension>::set_dof_values<BlockVector<double>,float>
(const Vector<float>&, BlockVector<double>&) const;
template
void
DoFObjectAccessor<2,deal_II_dimension>::get_dof_values<BlockVector<float>,float>
(const BlockVector<float>&, Vector<float>&) const;
template
void
DoFObjectAccessor<2,deal_II_dimension>::set_dof_values<BlockVector<float>,float>
(const Vector<float>&, BlockVector<float>&) const;


#endif


#if deal_II_dimension >= 3
// for double
template
void
DoFObjectAccessor<3,deal_II_dimension>::get_dof_values<Vector<double>,double>
(const Vector<double>&, Vector<double>&) const;
template
void
DoFObjectAccessor<3,deal_II_dimension>::set_dof_values<Vector<double>,double>
(const Vector<double>&, Vector<double>&) const;
template
void
DoFObjectAccessor<3,deal_II_dimension>::get_dof_values<Vector<double>,float>
(const Vector<double>&, Vector<float>&) const;
template
void
DoFObjectAccessor<3,deal_II_dimension>::set_dof_values<Vector<float>,double>
(const Vector<double>&, Vector<float>&) const;
template
void
DoFObjectAccessor<3,deal_II_dimension>::get_dof_values<Vector<float>,double>
(const Vector<float>&, Vector<double>&) const;
template
void
DoFObjectAccessor<3,deal_II_dimension>::set_dof_values<Vector<double>,float>
(const Vector<float>&, Vector<double>&) const;
template
void
DoFObjectAccessor<3,deal_II_dimension>::get_dof_values<Vector<float>,float>
(const Vector<float>&, Vector<float>&) const;
template
void
DoFObjectAccessor<3,deal_II_dimension>::set_dof_values<Vector<float>,float>
(const Vector<float>&, Vector<float>&) const;



// for block vector
template
void
DoFObjectAccessor<3,deal_II_dimension>::get_dof_values<BlockVector<double>,double>
(const BlockVector<double>&, Vector<double>&) const;
template
void
DoFObjectAccessor<3,deal_II_dimension>::set_dof_values<BlockVector<double>,double>
(const Vector<double>&, BlockVector<double>&) const;
template
void
DoFObjectAccessor<3,deal_II_dimension>::get_dof_values<BlockVector<double>,float>
(const BlockVector<double>&, Vector<float>&) const;
template
void
DoFObjectAccessor<3,deal_II_dimension>::set_dof_values<BlockVector<float>,double>
(const Vector<double>&, BlockVector<float>&) const;
template
void
DoFObjectAccessor<3,deal_II_dimension>::get_dof_values<BlockVector<float>,double>
(const BlockVector<float>&, Vector<double>&) const;
template
void
DoFObjectAccessor<3,deal_II_dimension>::set_dof_values<BlockVector<double>,float>
(const Vector<float>&, BlockVector<double>&) const;
template
void
DoFObjectAccessor<3,deal_II_dimension>::get_dof_values<BlockVector<float>,float>
(const BlockVector<float>&, Vector<float>&) const;
template
void
DoFObjectAccessor<3,deal_II_dimension>::set_dof_values<BlockVector<float>,float>
(const Vector<float>&, BlockVector<float>&) const;


#endif



template
void
DoFCellAccessor<deal_II_dimension>::
get_interpolated_dof_values<Vector<double>,double>
(const Vector<double>&, Vector<double>&) const;
template
void
DoFCellAccessor<deal_II_dimension>::
set_dof_values_by_interpolation<Vector<double>,double>
(const Vector<double>&, Vector<double>&) const;

template
void
DoFCellAccessor<deal_II_dimension>::
get_interpolated_dof_values<Vector<double>,float>
(const Vector<double>&, Vector<float>&) const;
template
void
DoFCellAccessor<deal_II_dimension>::
set_dof_values_by_interpolation<Vector<float>,double>
(const Vector<double>&, Vector<float>&) const;

template
void
DoFCellAccessor<deal_II_dimension>::
get_interpolated_dof_values<Vector<float>,double>
(const Vector<float>&, Vector<double>&) const;
template
void
DoFCellAccessor<deal_II_dimension>::
set_dof_values_by_interpolation<Vector<double>,float>
(const Vector<float>&, Vector<double>&) const;

template
void
DoFCellAccessor<deal_II_dimension>::
get_interpolated_dof_values<Vector<float>,float>
(const Vector<float>&, Vector<float>&) const;
template
void
DoFCellAccessor<deal_II_dimension>::
set_dof_values_by_interpolation<Vector<float>,float>
(const Vector<float>&, Vector<float>&) const;


template
void
DoFCellAccessor<deal_II_dimension>::
get_interpolated_dof_values<BlockVector<double>,double>
(const BlockVector<double>&, Vector<double>&) const;

template
void
DoFCellAccessor<deal_II_dimension>::
set_dof_values_by_interpolation<BlockVector<double>,double>
(const Vector<double>&, BlockVector<double>&) const;

template
void
DoFCellAccessor<deal_II_dimension>::
get_interpolated_dof_values<BlockVector<double>,float>
(const BlockVector<double>&, Vector<float>&) const;

template
void
DoFCellAccessor<deal_II_dimension>::
set_dof_values_by_interpolation<BlockVector<float>,double>
(const Vector<double>&, BlockVector<float>&) const;

template
void
DoFCellAccessor<deal_II_dimension>::
get_interpolated_dof_values<BlockVector<float>,double>
(const BlockVector<float>&, Vector<double>&) const;
template
void
DoFCellAccessor<deal_II_dimension>::
set_dof_values_by_interpolation<BlockVector<double>,float>
(const Vector<float>&, BlockVector<double>&) const;

template
void
DoFCellAccessor<deal_II_dimension>::
get_interpolated_dof_values<BlockVector<float>,float>
(const BlockVector<float>&, Vector<float>&) const;
template
void
DoFCellAccessor<deal_II_dimension>::
set_dof_values_by_interpolation<BlockVector<float>,float>
(const Vector<float>&, BlockVector<float>&) const;


template class DoFAccessor<deal_II_dimension>;

#if deal_II_dimension == 1
template class DoFObjectAccessor<1, 1>;
#endif

#if deal_II_dimension == 2
template class DoFObjectAccessor<1, 2>;
template class DoFObjectAccessor<2, 2>;

template class TriaRawIterator   <2,DoFObjectAccessor<1, 2> >;
template class TriaIterator      <2,DoFObjectAccessor<1, 2> >;
template class TriaActiveIterator<2,DoFObjectAccessor<1, 2> >;
#endif


#if deal_II_dimension == 3
template class DoFObjectAccessor<1, 3>;
template class DoFObjectAccessor<2, 3>;
template class DoFObjectAccessor<3, 3>;

template class TriaRawIterator   <3,DoFObjectAccessor<1, 3> >;
template class TriaIterator      <3,DoFObjectAccessor<1, 3> >;
template class TriaActiveIterator<3,DoFObjectAccessor<1, 3> >;
template class TriaRawIterator   <3,DoFObjectAccessor<2, 3> >;
template class TriaIterator      <3,DoFObjectAccessor<2, 3> >;
template class TriaActiveIterator<3,DoFObjectAccessor<2, 3> >;
#endif


template class DoFCellAccessor<deal_II_dimension>;

template class TriaRawIterator   <deal_II_dimension,DoFCellAccessor<deal_II_dimension> >;
template class TriaIterator      <deal_II_dimension,DoFCellAccessor<deal_II_dimension> >;
template class TriaActiveIterator<deal_II_dimension,DoFCellAccessor<deal_II_dimension> >;
