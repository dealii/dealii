//----------------------------  mg_dof_accessor.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  mg_dof_accessor.cc  ---------------------------


#include <dofs/dof_levels.h>
#include <multigrid/mg_dof_accessor.h>
#include <multigrid/mg_dof_handler.h>
#include <grid/tria_iterator.h>
#include <grid/tria_iterator.templates.h>
#include <fe/fe.h>

#include <lac/vector.h>
#include <lac/full_matrix.h>
#include <lac/sparse_matrix.h>



/* ------------------------ MGDoFLineAccessor --------------------------- */

template <int dim>
MGDoFObjectAccessor<1, dim>::MGDoFObjectAccessor (const Triangulation<dim> *tria,
						  const int                 level,
						  const int                 index,
						  const AccessorData       *local_data) :
		MGDoFAccessor<dim> (local_data),
		MGDoFObjectAccessor_Inheritance<1,dim>::BaseClass(tria,level,index,local_data)
{}


template <int dim>
unsigned int MGDoFObjectAccessor<1, dim>::mg_dof_index (const unsigned int i) const
{
  Assert (this->dof_handler != 0,
	  typename DoFAccessor<dim>::ExcInvalidObject());
  Assert (this->mg_dof_handler != 0,
	  typename DoFAccessor<dim>::ExcInvalidObject());
				   // make sure a FE has been selected
				   // and enough room was reserved
  Assert (&this->dof_handler->get_fe() != 0,
	  typename DoFAccessor<dim>::ExcInvalidObject());
  Assert (i<this->dof_handler->get_fe().dofs_per_line,
	  ExcIndexRange (i, 0, this->dof_handler->get_fe().dofs_per_line));

  return this->mg_dof_handler->mg_levels[this->present_level]
    ->line_dofs[this->present_index*this->dof_handler->get_fe().dofs_per_line+i];
}


template <int dim>
void MGDoFObjectAccessor<1, dim>::set_mg_dof_index (const unsigned int i,
						    const unsigned int index) const
{
  Assert (this->dof_handler != 0,
	  typename DoFAccessor<dim>::ExcInvalidObject());
  Assert (this->mg_dof_handler != 0,
	  typename DoFAccessor<dim>::ExcInvalidObject());
				   // make sure a FE has been selected
				   // and enough room was reserved
  Assert (&this->dof_handler->get_fe() != 0,
	  typename DoFAccessor<dim>::ExcInvalidObject());
  Assert (i<this->dof_handler->get_fe().dofs_per_line,
	  ExcIndexRange (i, 0, this->dof_handler->get_fe().dofs_per_line));

  this->mg_dof_handler->mg_levels[this->present_level]
    ->line_dofs[this->present_index*this->dof_handler->get_fe().dofs_per_line+i] = index;
}


template <int dim>
unsigned int MGDoFObjectAccessor<1, dim>::mg_vertex_dof_index (const unsigned int vertex,
							       const unsigned int i) const
{
  Assert (this->dof_handler != 0,
	  typename DoFAccessor<dim>::ExcInvalidObject());
  Assert (this->mg_dof_handler != 0,
	  typename DoFAccessor<dim>::ExcInvalidObject());
  Assert (&this->dof_handler->get_fe() != 0,
	  typename DoFAccessor<dim>::ExcInvalidObject());
  Assert (vertex<2, ExcIndexRange (i,0,2));
  Assert (i<this->dof_handler->get_fe().dofs_per_vertex,
	  ExcIndexRange (i, 0, this->dof_handler->get_fe().dofs_per_vertex));

  return (this->mg_dof_handler->mg_vertex_dofs[this->vertex_index(vertex)]
	  .get_index (this->present_level, i, this->dof_handler->get_fe().dofs_per_vertex));
}


template <int dim>
void MGDoFObjectAccessor<1, dim>::set_mg_vertex_dof_index (const unsigned int vertex,
							   const unsigned int i,
							   const unsigned int index) const
{
  Assert (this->dof_handler != 0,
	  typename DoFAccessor<dim>::ExcInvalidObject());
  Assert (this->mg_dof_handler != 0,
	  typename DoFAccessor<dim>::ExcInvalidObject());
  Assert (&this->dof_handler->get_fe() != 0,
	  typename DoFAccessor<dim>::ExcInvalidObject());
  Assert (vertex<2, ExcIndexRange (i,0,2));
  Assert (i<this->dof_handler->get_fe().dofs_per_vertex,
	  ExcIndexRange (i, 0, this->dof_handler->get_fe().dofs_per_vertex));

  this->mg_dof_handler->mg_vertex_dofs[this->vertex_index(vertex)]
    .set_index (this->present_level, i, this->dof_handler->get_fe().dofs_per_vertex, index);
}


template <int dim>
void
MGDoFObjectAccessor<1, dim>::get_mg_dof_indices (std::vector<unsigned int> &dof_indices) const
{
  Assert (this->dof_handler != 0,
	  typename DoFAccessor<dim>::ExcInvalidObject());
  Assert (this->mg_dof_handler != 0,
	  typename DoFAccessor<dim>::ExcInvalidObject());
  Assert (&this->dof_handler->get_fe() != 0,
	  typename DoFAccessor<dim>::ExcInvalidObject());
  Assert (dof_indices.size() == (2*this->dof_handler->get_fe().dofs_per_vertex +
				 this->dof_handler->get_fe().dofs_per_line),
	  typename DoFAccessor<dim>::ExcVectorDoesNotMatch());

  const unsigned int dofs_per_vertex = this->dof_handler->get_fe().dofs_per_vertex,
		     dofs_per_line   = this->dof_handler->get_fe().dofs_per_line;
  std::vector<unsigned int>::iterator next = dof_indices.begin();
  for (unsigned int vertex=0; vertex<2; ++vertex)
    for (unsigned int d=0; d<dofs_per_vertex; ++d)
      *next++ = mg_vertex_dof_index(vertex,d);
  for (unsigned int d=0; d<dofs_per_line; ++d)
    *next++ = mg_dof_index(d);
  
  Assert (next == dof_indices.end(),
	  ExcInternalError());
}


template <int dim>
template <typename number>
void
MGDoFObjectAccessor<1,dim>::get_mg_dof_values (const Vector<number> &values,
					       Vector<number>       &dof_values) const
{
  Assert (this->dof_handler != 0, DoFAccessor<1>::ExcInvalidObject());
  Assert (this->mg_dof_handler != 0, DoFAccessor<1>::ExcInvalidObject());
  Assert (&this->dof_handler->get_fe() != 0, DoFAccessor<1>::ExcInvalidObject());
  Assert (dof_values.size() == this->dof_handler->get_fe().dofs_per_cell,
	  DoFAccessor<1>::ExcVectorDoesNotMatch());
  Assert (values.size() == this->dof_handler->n_dofs(),
	  DoFAccessor<1>::ExcVectorDoesNotMatch());

  const unsigned int dofs_per_vertex = this->dof_handler->get_fe().dofs_per_vertex,
		     dofs_per_line   = this->dof_handler->get_fe().dofs_per_line;
  typename Vector<number>::iterator next_dof_value=dof_values.begin();
  for (unsigned int vertex=0; vertex<2; ++vertex)
    for (unsigned int d=0; d<dofs_per_vertex; ++d)
      *next_dof_value++ = values(mg_vertex_dof_index(vertex,d));
  for (unsigned int d=0; d<dofs_per_line; ++d)
    *next_dof_value++ = values(mg_dof_index(d));
  
  Assert (next_dof_value == dof_values.end(),
	  ExcInternalError());
}


template <int dim>
TriaIterator<dim,MGDoFObjectAccessor<1, dim> >
MGDoFObjectAccessor<1, dim>::child (const unsigned int i) const
{
  TriaIterator<dim,MGDoFObjectAccessor<1, dim> > q (this->tria,
						    this->present_level+1,
						    this->child_index (i),
						    this->mg_dof_handler);
  
#ifdef DEBUG
  if (q.state() != IteratorState::past_the_end)
    Assert (q->used(), typename TriaAccessor<dim>::ExcUnusedCellAsChild());
#endif
  return q;
}


template <int dim>
void
MGDoFObjectAccessor<1, dim>::copy_from (const MGDoFObjectAccessor<1, dim> &a)
{
  DoFObjectAccessor<1, dim>::copy_from (a);
  this->set_mg_dof_handler (a.mg_dof_handler);
}


/* ------------------------ MGDoFQuadAccessor --------------------------- */

template <int dim>
MGDoFObjectAccessor<2, dim>::MGDoFObjectAccessor (const Triangulation<dim> *tria,
						  const int                 level,
						  const int                 index,
						  const AccessorData       *local_data) :
		MGDoFAccessor<dim> (local_data),
		MGDoFObjectAccessor_Inheritance<2,dim>::BaseClass(tria,level,index,local_data)
{}


template <int dim>
unsigned int MGDoFObjectAccessor<2, dim>::mg_dof_index (const unsigned int i) const
{
  Assert (this->dof_handler != 0,
	  typename DoFAccessor<dim>::ExcInvalidObject());
  Assert (this->mg_dof_handler != 0,
	  typename DoFAccessor<dim>::ExcInvalidObject());
				   // make sure a FE has been selected
				   // and enough room was reserved
  Assert (&this->dof_handler->get_fe() != 0,
	  typename DoFAccessor<dim>::ExcInvalidObject());
  Assert (i<this->dof_handler->get_fe().dofs_per_quad,
	  ExcIndexRange (i, 0, this->dof_handler->get_fe().dofs_per_quad));

  return this->mg_dof_handler->mg_levels[this->present_level]
    ->quad_dofs[this->present_index*this->dof_handler->get_fe().dofs_per_quad+i];
}


template <int dim>
void MGDoFObjectAccessor<2, dim>::set_mg_dof_index (const unsigned int i,
						    const unsigned int index) const
{
  Assert (this->dof_handler != 0,
	  typename DoFAccessor<dim>::ExcInvalidObject());
  Assert (this->mg_dof_handler != 0,
	  typename DoFAccessor<dim>::ExcInvalidObject());
				   // make sure a FE has been selected
				   // and enough room was reserved
  Assert (&this->dof_handler->get_fe() != 0,
	  typename DoFAccessor<dim>::ExcInvalidObject());
  Assert (i<this->dof_handler->get_fe().dofs_per_quad,
	  ExcIndexRange (i, 0, this->dof_handler->get_fe().dofs_per_quad));

  this->mg_dof_handler->mg_levels[this->present_level]
    ->quad_dofs[this->present_index*this->dof_handler->get_fe().dofs_per_quad+i] = index;
}


template <int dim>
unsigned int MGDoFObjectAccessor<2, dim>::mg_vertex_dof_index (const unsigned int vertex,
							       const unsigned int i) const
{
  Assert (this->dof_handler != 0,
	  typename DoFAccessor<dim>::ExcInvalidObject());
  Assert (this->mg_dof_handler != 0,
	  typename DoFAccessor<dim>::ExcInvalidObject());
  Assert (&this->dof_handler->get_fe() != 0,
	  typename DoFAccessor<dim>::ExcInvalidObject());
  Assert (vertex<4, ExcIndexRange (i,0,4));
  Assert (i<this->dof_handler->get_fe().dofs_per_vertex,
	  ExcIndexRange (i, 0, this->dof_handler->get_fe().dofs_per_vertex));
  
  return (this->mg_dof_handler->mg_vertex_dofs[this->vertex_index(vertex)]
	  .get_index (this->present_level, i, this->dof_handler->get_fe().dofs_per_vertex));
}


template <int dim>
void MGDoFObjectAccessor<2, dim>::set_mg_vertex_dof_index (const unsigned int vertex,
							   const unsigned int i,
							   const unsigned int index) const
{
  Assert (this->dof_handler != 0,
	  typename DoFAccessor<dim>::ExcInvalidObject());
  Assert (this->mg_dof_handler != 0,
	  typename DoFAccessor<dim>::ExcInvalidObject());
  Assert (&this->dof_handler->get_fe() != 0,
	  typename DoFAccessor<dim>::ExcInvalidObject());
  Assert (vertex<4, ExcIndexRange (i,0,4));
  Assert (i<this->dof_handler->get_fe().dofs_per_vertex,
	  ExcIndexRange (i, 0, this->dof_handler->get_fe().dofs_per_vertex));

  this->mg_dof_handler->mg_vertex_dofs[this->vertex_index(vertex)]
    .set_index (this->present_level, i, this->dof_handler->get_fe().dofs_per_vertex, index);
}


template <int dim>
void
MGDoFObjectAccessor<2, dim>::get_mg_dof_indices (std::vector<unsigned int> &dof_indices) const
{
  Assert (this->dof_handler != 0,
	  typename DoFAccessor<dim>::ExcInvalidObject());
  Assert (this->mg_dof_handler != 0,
	  typename DoFAccessor<dim>::ExcInvalidObject());
  Assert (&this->dof_handler->get_fe() != 0,
	  typename DoFAccessor<dim>::ExcInvalidObject());
  Assert (dof_indices.size() == (4*this->dof_handler->get_fe().dofs_per_vertex +
				 4*this->dof_handler->get_fe().dofs_per_line +
				 this->dof_handler->get_fe().dofs_per_quad),
	  DoFAccessor<2>::ExcVectorDoesNotMatch());

  const unsigned int dofs_per_vertex = this->dof_handler->get_fe().dofs_per_vertex,
		     dofs_per_line   = this->dof_handler->get_fe().dofs_per_line,
		     dofs_per_quad   = this->dof_handler->get_fe().dofs_per_quad;
  std::vector<unsigned int>::iterator next = dof_indices.begin();
  for (unsigned int vertex=0; vertex<4; ++vertex)
    for (unsigned int d=0; d<dofs_per_vertex; ++d)
      *next++ = mg_vertex_dof_index(vertex,d);
  for (unsigned int line=0; line<4; ++line)
    for (unsigned int d=0; d<dofs_per_line; ++d)
      *next++ = this->line(line)->mg_dof_index(d);
  for (unsigned int d=0; d<dofs_per_quad; ++d)
    *next++ = mg_dof_index(d);
  
  Assert (next == dof_indices.end(),
	  ExcInternalError());
}


template <int dim>
template <typename number>
void
MGDoFObjectAccessor<2,dim>::get_mg_dof_values (const Vector<number> &values,
					       Vector<number>       &dof_values) const
{
  Assert (this->dof_handler != 0, DoFAccessor<2>::ExcInvalidObject());
  Assert (this->mg_dof_handler != 0, DoFAccessor<2>::ExcInvalidObject());
  Assert (&this->dof_handler->get_fe() != 0, DoFAccessor<2>::ExcInvalidObject());
  Assert (dof_values.size() == this->dof_handler->get_fe().dofs_per_cell,
	  DoFAccessor<2>::ExcVectorDoesNotMatch());
  Assert (values.size() == this->mg_dof_handler->n_dofs(this->present_level),
	  DoFAccessor<2>::ExcVectorDoesNotMatch());

  const unsigned int dofs_per_vertex = this->dof_handler->get_fe().dofs_per_vertex,
		     dofs_per_line   = this->dof_handler->get_fe().dofs_per_line,
		     dofs_per_quad   = this->dof_handler->get_fe().dofs_per_quad;
  typename Vector<number>::iterator next_dof_value=dof_values.begin();
  for (unsigned int vertex=0; vertex<4; ++vertex)
    for (unsigned int d=0; d<dofs_per_vertex; ++d)
      *next_dof_value++ = values(mg_vertex_dof_index(vertex,d));
  for (unsigned int line=0; line<4; ++line)
    for (unsigned int d=0; d<dofs_per_line; ++d)
      *next_dof_value++ = values(this->line(line)->mg_dof_index(d));
  for (unsigned int d=0; d<dofs_per_quad; ++d)
    *next_dof_value++ = values(mg_dof_index(d));
  
  Assert (next_dof_value == dof_values.end(),
	  ExcInternalError());
}


template <int dim>
TriaIterator<dim,MGDoFObjectAccessor<1, dim> >
MGDoFObjectAccessor<2, dim>::line (const unsigned int i) const
{
  Assert (i<4, ExcIndexRange (i, 0, 4));

  return TriaIterator<dim,MGDoFObjectAccessor<1, dim> >
    (
      this->tria,
      this->present_level,
      this->line_index (i),
      this->mg_dof_handler
    );
}


template <int dim>
TriaIterator<dim,MGDoFObjectAccessor<2, dim> >
MGDoFObjectAccessor<2, dim>::child (const unsigned int i) const
{
  TriaIterator<dim,MGDoFObjectAccessor<2, dim> > q (this->tria,
						    this->present_level+1,
						    this->child_index (i),
						    this->mg_dof_handler);
  
#ifdef DEBUG
  if (q.state() != IteratorState::past_the_end)
    Assert (q->used(), typename TriaAccessor<dim>::ExcUnusedCellAsChild());
#endif
  return q;
}


template <int dim>
void
MGDoFObjectAccessor<2, dim>::copy_from (const MGDoFObjectAccessor<2, dim> &a)
{
  DoFObjectAccessor<2, dim>::copy_from (a);
  this->set_mg_dof_handler (a.mg_dof_handler);
}


/* ------------------------ MGDoFHexAccessor --------------------------- */

template <int dim>
MGDoFObjectAccessor<3, dim>::MGDoFObjectAccessor (const Triangulation<dim> *tria,
						  const int                 level,
						  const int                 index,
						  const AccessorData       *local_data) :
		MGDoFAccessor<dim> (local_data),
		MGDoFObjectAccessor_Inheritance<3,dim>::BaseClass(tria,level,index,local_data)
{}


template <int dim>
unsigned int MGDoFObjectAccessor<3, dim>::mg_dof_index (const unsigned int i) const
{
  Assert (this->dof_handler != 0,
	  typename DoFAccessor<dim>::ExcInvalidObject());
  Assert (this->mg_dof_handler != 0,
	  typename DoFAccessor<dim>::ExcInvalidObject());
				   // make sure a FE has been selected
				   // and enough room was reserved
  Assert (&this->dof_handler->get_fe() != 0,
	  typename DoFAccessor<dim>::ExcInvalidObject());
  Assert (i<this->dof_handler->get_fe().dofs_per_hex,
	  ExcIndexRange (i, 0, this->dof_handler->get_fe().dofs_per_hex));

  return this->mg_dof_handler->mg_levels[this->present_level]
    ->hex_dofs[this->present_index*this->dof_handler->get_fe().dofs_per_hex+i];
}


template <int dim>
void MGDoFObjectAccessor<3, dim>::set_mg_dof_index (const unsigned int i,
						    const unsigned int index) const
{
  Assert (this->dof_handler != 0,
	  typename DoFAccessor<dim>::ExcInvalidObject());
  Assert (this->mg_dof_handler != 0,
	  typename DoFAccessor<dim>::ExcInvalidObject());
				   // make sure a FE has been selected
				   // and enough room was reserved
  Assert (&this->dof_handler->get_fe() != 0,
	  typename DoFAccessor<dim>::ExcInvalidObject());
  Assert (i<this->dof_handler->get_fe().dofs_per_hex,
	  ExcIndexRange (i, 0, this->dof_handler->get_fe().dofs_per_hex));

  this->mg_dof_handler->mg_levels[this->present_level]
    ->hex_dofs[this->present_index*this->dof_handler->get_fe().dofs_per_hex+i] = index;
}


template <int dim>
unsigned int MGDoFObjectAccessor<3, dim>::mg_vertex_dof_index (const unsigned int vertex,
							       const unsigned int i) const
{
  Assert (this->dof_handler != 0,
	  typename DoFAccessor<dim>::ExcInvalidObject());
  Assert (this->mg_dof_handler != 0,
	  typename DoFAccessor<dim>::ExcInvalidObject());
  Assert (&this->dof_handler->get_fe() != 0,
	  typename DoFAccessor<dim>::ExcInvalidObject());
  Assert (vertex<8, ExcIndexRange (i,0,8));
  Assert (i<this->dof_handler->get_fe().dofs_per_vertex,
	  ExcIndexRange (i, 0, this->dof_handler->get_fe().dofs_per_vertex));
  
  return (this->mg_dof_handler->mg_vertex_dofs[this->vertex_index(vertex)]
	  .get_index (this->present_level, i, this->dof_handler->get_fe().dofs_per_vertex));
}


template <int dim>
void MGDoFObjectAccessor<3, dim>::set_mg_vertex_dof_index (const unsigned int vertex,
							   const unsigned int i,
							   const unsigned int index) const
{
  Assert (this->dof_handler != 0,
	  typename DoFAccessor<dim>::ExcInvalidObject());
  Assert (this->mg_dof_handler != 0,
	  typename DoFAccessor<dim>::ExcInvalidObject());
  Assert (&this->dof_handler->get_fe() != 0,
	  typename DoFAccessor<dim>::ExcInvalidObject());
  Assert (vertex<4, ExcIndexRange (i,0,4));
  Assert (i<this->dof_handler->get_fe().dofs_per_vertex,
	  ExcIndexRange (i, 0, this->dof_handler->get_fe().dofs_per_vertex));

  this->mg_dof_handler->mg_vertex_dofs[this->vertex_index(vertex)]
    .set_index (this->present_level, i, this->dof_handler->get_fe().dofs_per_vertex, index);
}


template <int dim>
void
MGDoFObjectAccessor<3, dim>::get_mg_dof_indices (std::vector<unsigned int> &dof_indices) const
{
  Assert (this->dof_handler != 0,
	  typename DoFAccessor<dim>::ExcInvalidObject());
  Assert (this->mg_dof_handler != 0,
	  typename DoFAccessor<dim>::ExcInvalidObject());
  Assert (&this->dof_handler->get_fe() != 0,
	  typename DoFAccessor<dim>::ExcInvalidObject());
  Assert (dof_indices.size() == (8*this->dof_handler->get_fe().dofs_per_vertex +
				 12*this->dof_handler->get_fe().dofs_per_line +
				 6*this->dof_handler->get_fe().dofs_per_quad +
				 this->dof_handler->get_fe().dofs_per_hex),
	  DoFAccessor<3>::ExcVectorDoesNotMatch());

  const unsigned int dofs_per_vertex = this->dof_handler->get_fe().dofs_per_vertex,
		     dofs_per_line   = this->dof_handler->get_fe().dofs_per_line,
		     dofs_per_quad   = this->dof_handler->get_fe().dofs_per_quad,
		     dofs_per_hex    = this->dof_handler->get_fe().dofs_per_hex;
  std::vector<unsigned int>::iterator next = dof_indices.begin();
  for (unsigned int vertex=0; vertex<8; ++vertex)
    for (unsigned int d=0; d<dofs_per_vertex; ++d)
      *next++ = mg_vertex_dof_index(vertex,d);
  for (unsigned int line=0; line<12; ++line)
    for (unsigned int d=0; d<dofs_per_line; ++d)
      *next++ = this->line(line)->mg_dof_index(d);
  for (unsigned int quad=0; quad<12; ++quad)
    for (unsigned int d=0; d<dofs_per_quad; ++d)
      *next++ = this->quad(quad)->mg_dof_index(d);
  for (unsigned int d=0; d<dofs_per_hex; ++d)
    *next++ = mg_dof_index(d);
  
  Assert (next == dof_indices.end(),
	  ExcInternalError());
}


template <int dim>
template <typename number>
void
MGDoFObjectAccessor<3,dim>::get_mg_dof_values (const Vector<number> &values,
					       Vector<number>       &dof_values) const
{
  Assert (this->dof_handler != 0, DoFAccessor<2>::ExcInvalidObject());
  Assert (this->mg_dof_handler != 0, DoFAccessor<2>::ExcInvalidObject());
  Assert (&this->dof_handler->get_fe() != 0, DoFAccessor<2>::ExcInvalidObject());
  Assert (dof_values.size() == this->dof_handler->get_fe().dofs_per_cell,
	  DoFAccessor<3>::ExcVectorDoesNotMatch());
  Assert (values.size() == this->mg_dof_handler->n_dofs(this->present_level),
	  DoFAccessor<3>::ExcVectorDoesNotMatch());

  const unsigned int dofs_per_vertex = this->dof_handler->get_fe().dofs_per_vertex,
		     dofs_per_line   = this->dof_handler->get_fe().dofs_per_line,
		     dofs_per_quad   = this->dof_handler->get_fe().dofs_per_quad,
		     dofs_per_hex    = this->dof_handler->get_fe().dofs_per_hex;
  typename Vector<number>::iterator next_dof_value=dof_values.begin();
  for (unsigned int vertex=0; vertex<8; ++vertex)
    for (unsigned int d=0; d<dofs_per_vertex; ++d)
      *next_dof_value++ = values(mg_vertex_dof_index(vertex,d));
  for (unsigned int line=0; line<12; ++line)
    for (unsigned int d=0; d<dofs_per_line; ++d)
      *next_dof_value++ = values(this->line(line)->mg_dof_index(d));
  for (unsigned int quad=0; quad<6; ++quad)
    for (unsigned int d=0; d<dofs_per_quad; ++d)
      *next_dof_value++ = values(this->quad(quad)->mg_dof_index(d));
  for (unsigned int d=0; d<dofs_per_hex; ++d)
    *next_dof_value++ = values(this->dof_index(d));
  
  Assert (next_dof_value == dof_values.end(),
	  ExcInternalError());
}


template <int dim>
TriaIterator<dim,MGDoFObjectAccessor<1, dim> >
MGDoFObjectAccessor<3, dim>::line (const unsigned int i) const {
  Assert (i<12, ExcIndexRange (i, 0, 12));

  return TriaIterator<dim,MGDoFObjectAccessor<1, dim> >
    (
      this->tria,
      this->present_level,
      this->line_index (i),
      this->mg_dof_handler
    );
}


template <int dim>
TriaIterator<dim,MGDoFObjectAccessor<2, dim> >
MGDoFObjectAccessor<3, dim>::quad (const unsigned int i) const {
  Assert (i<12, ExcIndexRange (i, 0, 6));

  return TriaIterator<dim,MGDoFObjectAccessor<2, dim> >
    (
      this->tria,
      this->present_level,
      this->quad_index (i),
      this->mg_dof_handler
    );
}


template <int dim>
TriaIterator<dim,MGDoFObjectAccessor<3, dim> >
MGDoFObjectAccessor<3, dim>::child (const unsigned int i) const {
  TriaIterator<dim,MGDoFObjectAccessor<3, dim> > q (this->tria,
						    this->present_level+1,
						    this->child_index (i),
						    this->mg_dof_handler);
  
#ifdef DEBUG
  if (q.state() != IteratorState::past_the_end)
    Assert (q->used(), typename TriaAccessor<dim>::ExcUnusedCellAsChild());
#endif
  return q;
}


template <int dim>
void
MGDoFObjectAccessor<3, dim>::copy_from (const MGDoFObjectAccessor<3, dim> &a) {
  DoFObjectAccessor<3, dim>::copy_from (a);
  this->set_mg_dof_handler (a.mg_dof_handler);
}


/*------------------------- Functions: MGDoFCellAccessor -----------------------*/

template <int dim>
TriaIterator<dim,MGDoFCellAccessor<dim> >
MGDoFCellAccessor<dim>::neighbor (const unsigned int i) const {
  TriaIterator<dim,MGDoFCellAccessor<dim> > q (this->tria,
					       this->neighbor_level (i),
					       this->neighbor_index (i),
					       this->mg_dof_handler);
  
#ifdef DEBUG
  if (q.state() != IteratorState::past_the_end)
    Assert (q->used(), typename TriaAccessor<dim>::ExcUnusedCellAsNeighbor());
#endif
  return q;
}


template <int dim>
TriaIterator<dim,MGDoFCellAccessor<dim> >
MGDoFCellAccessor<dim>::child (const unsigned int i) const {
  TriaIterator<dim,MGDoFCellAccessor<dim> > q (this->tria,
					       this->present_level+1,
					       this->child_index (i),
					       this->mg_dof_handler);
  
#ifdef DEBUG
  if (q.state() != IteratorState::past_the_end)
    Assert (q->used(), typename TriaAccessor<dim>::ExcUnusedCellAsChild());
#endif
  return q;
}


#if deal_II_dimension == 1

template <>
MGDoFCellAccessor<1>::face_iterator
MGDoFCellAccessor<1>::face (const unsigned int) const
{
  Assert (false, ExcNotUsefulForThisDimension());
  return face_iterator();
}

#endif


#if deal_II_dimension == 2

template <>
MGDoFCellAccessor<2>::face_iterator
MGDoFCellAccessor<2>::face (const unsigned int i) const
{
  return this->line(i);
}

#endif


#if deal_II_dimension == 3

template <>
MGDoFCellAccessor<3>::face_iterator
MGDoFCellAccessor<3>::face (const unsigned int i) const
{
  return this->quad(i);
}

#endif


// explicit instantiations

template
void
MGDoFObjectAccessor<1,deal_II_dimension>::
get_mg_dof_values (const Vector<double> &values,
		   Vector<double>       &dof_values) const;

template
void
MGDoFObjectAccessor<1,deal_II_dimension>::
get_mg_dof_values (const Vector<float> &values,
		   Vector<float>       &dof_values) const;


#if deal_II_dimension >= 2

template
void
MGDoFObjectAccessor<2,deal_II_dimension>::
get_mg_dof_values (const Vector<double> &values,
		   Vector<double>       &dof_values) const;

template
void
MGDoFObjectAccessor<2,deal_II_dimension>::
get_mg_dof_values (const Vector<float> &values,
		   Vector<float>       &dof_values) const;

#endif


#if deal_II_dimension >= 3

template
void
MGDoFObjectAccessor<3,deal_II_dimension>::
get_mg_dof_values (const Vector<double> &values,
		   Vector<double>       &dof_values) const;

template
void
MGDoFObjectAccessor<3,deal_II_dimension>::
get_mg_dof_values (const Vector<float> &values,
		   Vector<float>       &dof_values) const;

#endif


#if deal_II_dimension == 1
template class MGDoFObjectAccessor<1, 1>;
template class MGDoFCellAccessor<1>;

template class TriaRawIterator<1,MGDoFCellAccessor<1> >;
template class TriaIterator<1,MGDoFCellAccessor<1> >;
template class TriaActiveIterator<1,MGDoFCellAccessor<1> >;
#endif

#if deal_II_dimension == 2
template class MGDoFObjectAccessor<1, 2>;
template class MGDoFObjectAccessor<2, 2>;
template class MGDoFCellAccessor<2>;

template class TriaRawIterator   <2,MGDoFObjectAccessor<1, 2> >;
template class TriaIterator      <2,MGDoFObjectAccessor<1, 2> >;
template class TriaActiveIterator<2,MGDoFObjectAccessor<1, 2> >;

template class TriaRawIterator   <2,MGDoFCellAccessor<2> >;
template class TriaIterator      <2,MGDoFCellAccessor<2> >;
template class TriaActiveIterator<2,MGDoFCellAccessor<2> >;
#endif


#if deal_II_dimension == 3
template class MGDoFObjectAccessor<1, 3>;
template class MGDoFObjectAccessor<2, 3>;
template class MGDoFObjectAccessor<3, 3>;
template class MGDoFCellAccessor<3>;

template class TriaRawIterator   <3,MGDoFObjectAccessor<1, 3> >;
template class TriaIterator      <3,MGDoFObjectAccessor<1, 3> >;
template class TriaActiveIterator<3,MGDoFObjectAccessor<1, 3> >;

template class TriaRawIterator   <3,MGDoFObjectAccessor<2, 3> >;
template class TriaIterator      <3,MGDoFObjectAccessor<2, 3> >;
template class TriaActiveIterator<3,MGDoFObjectAccessor<2, 3> >;

template class TriaRawIterator   <3,MGDoFCellAccessor<3> >;
template class TriaIterator      <3,MGDoFCellAccessor<3> >;
template class TriaActiveIterator<3,MGDoFCellAccessor<3> >;
#endif
