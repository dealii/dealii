//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__dof_accessor_templates_h
#define __deal2__dof_accessor_templates_h


#include <base/config.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_levels.h>
#include <dofs/hp_dof_levels.h>
#include <grid/tria_iterator.h>
#include <grid/tria_iterator.templates.h>

#include <vector>


  

/*------------------------- Functions: DoFAccessor ---------------------------*/


template <class DH>
DoFAccessor<DH>::DoFAccessor ()
                :
		dof_handler(0)
{
  Assert (false, ExcInvalidObject());
}



template <class DH>
inline
DoFAccessor<DH>::DoFAccessor (const DH *dof_handler)
		:
		dof_handler(const_cast<DH*>(dof_handler))
{}



template <class DH>
void
DoFAccessor<DH>::set_dof_handler (DH *dh)
{
  Assert (dh != 0, ExcInvalidObject());
  this->dof_handler = dh;
}



template <class DH>
inline
const DH &
DoFAccessor<DH>::get_dof_handler () const
{
  return *this->dof_handler;
}


template <class DH>
inline
DoFAccessor<DH> &
DoFAccessor<DH>::operator = (const DoFAccessor<DH> &da)
{
  this->set_dof_handler (da.dof_handler);
  return *this;
}



template <class DH>
inline
bool
DoFAccessor<DH>::operator == (const DoFAccessor<DH> &a) const
{
  Assert (this->dof_handler == a.dof_handler, ExcCantCompareIterators());

                                   // there is no real data to compare, except
                                   // to make sure that the this->dof_handler
                                   // objects in use are the same
  return true;
}



template <class DH>
inline
bool
DoFAccessor<DH>::operator != (const DoFAccessor<DH> &a) const
{
  Assert (this->dof_handler == a.dof_handler, ExcCantCompareIterators());

                                   // there is no real data to compare, except
                                   // to make sure that the this->dof_handler
                                   // objects in use are the same. this is
                                   // checked above, and apart from this there
                                   // is no reason for us to believe that the
                                   // two accessors are different
  return false;
}




/*------------------------- Functions: DoFObjectAccessor<1,dim> -----------------------*/


template <class DH>
inline
DoFObjectAccessor<1,DH>::
DoFObjectAccessor (const Triangulation<DH::dimension> *tria,
                   const int                 level,
                   const int                 index,
                   const AccessorData       *local_data)
                :
                DoFAccessor<DH> (local_data),
                DoFObjectAccessor_Inheritance<1,DH::dimension>::BaseClass (tria,
									   level,
									   index)
{}



template <class DH>
inline
TriaIterator<DH::dimension,DoFObjectAccessor<1,DH> >
DoFObjectAccessor<1,DH>::child (const unsigned int i) const
{
  Assert (static_cast<unsigned int>(this->present_level) < this->dof_handler->levels.size(),
          ExcMessage ("DoFHandler not initialized"));

  TriaIterator<dim,DoFObjectAccessor<1,DH> > q (this->tria,
						this->present_level+1,
						this->child_index (i),
						this->dof_handler);
  
#ifdef DEBUG
  if (q.state() != IteratorState::past_the_end)
    Assert (q->used(), typename TriaAccessor<dim>::ExcUnusedCellAsChild());
#endif
  return q;
}



template <class DH>
inline
void
DoFObjectAccessor<1,DH>::copy_from (const DoFObjectAccessor<1,DH> &a)
{
  BaseClass::copy_from (a);
  this->set_dof_handler (a.dof_handler);
}



template <class DH>
inline
bool
DoFObjectAccessor<1,DH>::operator == (const DoFObjectAccessor<1,DH> &a) const
{
  return (TriaObjectAccessor<1,dim>::operator == (a)
          &&
          DoFAccessor<DH>::operator == (a));
}


template <class DH>
inline
bool
DoFObjectAccessor<1,DH>::operator != (const DoFObjectAccessor<1,DH> &a) const
{
  return (TriaObjectAccessor<1,dim>::operator != (a)
          ||
          DoFAccessor<DH>::operator != (a));
}


/*------------------------- Functions: DoFObjectAccessor<2,dim> -----------------------*/

template <class DH>
DoFObjectAccessor<2,DH>::
DoFObjectAccessor (const Triangulation<DH::dimension> *tria,
                   const int                 level,
                   const int                 index,
                   const AccessorData       *local_data)
                :
                DoFAccessor<DH> (local_data),
                DoFObjectAccessor_Inheritance<2,dim>::BaseClass (tria,
                                                                 level,
                                                                 index)
{}



template <class DH>
inline
TriaIterator<DH::dimension,DoFObjectAccessor<1,DH> >
DoFObjectAccessor<2,DH>::line (const unsigned int i) const
{
  Assert (i<4, ExcIndexRange (i, 0, 4));
  Assert (static_cast<unsigned int>(this->present_level) < this->dof_handler->levels.size(),
          ExcMessage ("DoFHandler not initialized"));

  return TriaIterator<dim,DoFObjectAccessor<1,DH> >
    (
      this->tria,
      this->present_level,
      this->line_index (i),
      this->dof_handler
      );
}



template <class DH>
inline
TriaIterator<DH::dimension,DoFObjectAccessor<2,DH> >
DoFObjectAccessor<2,DH>::child (const unsigned int i) const
{
  Assert (static_cast<unsigned int>(this->present_level) < this->dof_handler->levels.size(),
          ExcMessage ("DoFHandler not initialized"));

  TriaIterator<dim,DoFObjectAccessor<2,DH> > q (this->tria,
						this->present_level+1,
						this->child_index (i),
						this->dof_handler);
  
#ifdef DEBUG
  if (q.state() != IteratorState::past_the_end)
    Assert (q->used(), typename TriaAccessor<dim>::ExcUnusedCellAsChild());
#endif
  return q;
}



template <class DH>
inline
void
DoFObjectAccessor<2,DH>::copy_from (const DoFObjectAccessor<2,DH> &a)
{
  BaseClass::copy_from (a);
  this->set_dof_handler (a.dof_handler);
}



template <class DH>
inline
bool
DoFObjectAccessor<2,DH>::operator == (const DoFObjectAccessor<2,DH> &a) const
{
  return (TriaObjectAccessor<2,dim>::operator == (a)
          &&
          DoFAccessor<DH>::operator == (a));
}


template <class DH>
inline
bool
DoFObjectAccessor<2,DH>::operator != (const DoFObjectAccessor<2,DH> &a) const
{
  return (TriaObjectAccessor<2,dim>::operator != (a)
          ||
          DoFAccessor<DH>::operator != (a));
}



/*------------------------- Functions: DoFObjectAccessor<3,dim> -----------------------*/

template <class DH>
inline
DoFObjectAccessor<3,DH>::
DoFObjectAccessor (const Triangulation<DH::dimension> *tria,
                   const int                 level,
                   const int                 index,
                   const AccessorData       *local_data)
                :
                DoFAccessor<DH> (local_data),
                DoFObjectAccessor_Inheritance<3,dim>::BaseClass (tria,
                                                                 level,
                                                                 index)
{}



template <class DH>
inline
TriaIterator<DH::dimension,DoFObjectAccessor<1,DH> >
DoFObjectAccessor<3,DH>::line (const unsigned int i) const
{
  TriaIterator<dim,TriaObjectAccessor<1,dim> > l = BaseClass::line(i);
  return TriaIterator<dim,DoFObjectAccessor<1,DH> >
    (
      this->tria,
      this->present_level,
      l->index(),
      this->dof_handler
      );
}


template <class DH>
inline
TriaIterator<DH::dimension,DoFObjectAccessor<2,DH> >
DoFObjectAccessor<3,DH>::quad (const unsigned int i) const
{
  Assert (i<6, ExcIndexRange (i, 0, 6));

  return TriaIterator<dim,DoFObjectAccessor<2,DH> >
    (
      this->tria,
      this->present_level,
      this->quad_index (i),
      this->dof_handler
      );
}


template <class DH>
inline
TriaIterator<DH::dimension,DoFObjectAccessor<3,DH> >
DoFObjectAccessor<3,DH>::child (const unsigned int i) const
{
  TriaIterator<dim,DoFObjectAccessor<3,DH> > q (this->tria,
						this->present_level+1,
						this->child_index (i),
						this->dof_handler);
  
#ifdef DEBUG
  if (q.state() != IteratorState::past_the_end)
    Assert (q->used(), typename TriaAccessor<dim>::ExcUnusedCellAsChild());
#endif
  return q;
}


template <class DH>
void DoFObjectAccessor<3,DH>::copy_from (const DoFObjectAccessor<3,DH> &a)
{
  BaseClass::copy_from (a);
  this->set_dof_handler (a.dof_handler);
}



template <class DH>
inline
bool
DoFObjectAccessor<3,DH>::operator == (const DoFObjectAccessor<3,DH> &a) const
{
  return (TriaObjectAccessor<3,dim>::operator == (a)
          &&
          DoFAccessor<DH>::operator == (a));
}


template <class DH>
inline
bool
DoFObjectAccessor<3,DH>::operator != (const DoFObjectAccessor<3,DH> &a) const
{
  return (TriaObjectAccessor<3,dim>::operator != (a)
          ||
          DoFAccessor<DH>::operator != (a));
}


/*--------------- Functions: DoFObjectAccessor<1,DoFHandler> -----------*/


template <>
inline
unsigned int
DoFObjectAccessor<1,DoFHandler<1> >::dof_index (const unsigned int i,
						const unsigned int fe_index) const
{
  Assert (fe_index == DoFHandler<1>::default_fe_index,
	  ExcMessage ("Only the default FE index is allowed for non-hp DoFHandler objects"));
  Assert (static_cast<unsigned int>(this->present_level) < this->dof_handler->levels.size(),
          ExcMessage ("DoFHandler not initialized"));

  Assert (this->dof_handler != 0, DoFAccessor<DoFHandler<1> >::ExcInvalidObject());
				   // make sure a FE has been selected
				   // and enough room was reserved
  Assert (&this->dof_handler->get_fe() != 0, DoFAccessor<DoFHandler<1> >::ExcInvalidObject());
  Assert (i<this->dof_handler->get_fe().dofs_per_line,
	  ExcIndexRange (i, 0, this->dof_handler->get_fe().dofs_per_line));

  return this->dof_handler->levels[this->present_level]
    ->line_dofs[this->present_index*this->dof_handler->get_fe().dofs_per_line+i];
}



template <>
inline
unsigned int
DoFObjectAccessor<1,DoFHandler<1> >::vertex_dof_index (const unsigned int vertex,
						       const unsigned int i,
						       const unsigned int fe_index) const
{
  Assert (fe_index == DoFHandler<1>::default_fe_index,
	  ExcMessage ("Only the default FE index is allowed for non-hp DoFHandler objects"));
  Assert (static_cast<unsigned int>(this->present_level) < this->dof_handler->levels.size(),
          ExcMessage ("DoFHandler not initialized"));

  Assert (this->dof_handler != 0, DoFAccessor<DoFHandler<1> >::ExcInvalidObject());
  Assert (&this->dof_handler->get_fe() != 0, DoFAccessor<DoFHandler<1> >::ExcInvalidObject());
  Assert (vertex<2, ExcIndexRange (i,0,2));
  Assert (i<this->dof_handler->get_fe().dofs_per_vertex,
	  ExcIndexRange (i, 0, this->dof_handler->get_fe().dofs_per_vertex));

  const unsigned int dof_number = (this->vertex_index(vertex) *
				   this->dof_handler->get_fe().dofs_per_vertex +
				   i);
  return this->dof_handler->vertex_dofs[dof_number];
}



template <>
inline
void
DoFObjectAccessor<1,DoFHandler<1> >::get_dof_indices (std::vector<unsigned int> &dof_indices) const
{
  Assert (static_cast<unsigned int>(this->present_level) < this->dof_handler->levels.size(),
          ExcMessage ("DoFHandler not initialized"));

  Assert (this->dof_handler != 0, DoFAccessor<DoFHandler<1> >::ExcInvalidObject());
  Assert (&this->dof_handler->get_fe() != 0, DoFAccessor<DoFHandler<1> >::ExcInvalidObject());
  Assert (dof_indices.size() == (2*this->dof_handler->get_fe().dofs_per_vertex +
				 this->dof_handler->get_fe().dofs_per_line),
	  DoFAccessor<DoFHandler<1> >::ExcVectorDoesNotMatch());

				   // this function really only makes
				   // sense on non-active objects if
				   // all degrees of freedom are
				   // located on vertices, since
				   // otherwise there are degrees of
				   // freedom on sub-objects which are
				   // not allocated for this
				   // non-active thing
  Assert (!this->has_children() ||
	  (this->dof_handler->get_fe().dofs_per_cell ==
	   2*this->dof_handler->get_fe().dofs_per_vertex),
	  DoFAccessor<DoFHandler<1> >::ExcNotActive());
  
  const unsigned int dofs_per_vertex = this->dof_handler->get_fe().dofs_per_vertex,
		     dofs_per_line   = this->dof_handler->get_fe().dofs_per_line;
  std::vector<unsigned int>::iterator next = dof_indices.begin();
  for (unsigned int vertex=0; vertex<2; ++vertex)
    for (unsigned int d=0; d<dofs_per_vertex; ++d)
      *next++ = vertex_dof_index(vertex,d);
  for (unsigned int d=0; d<dofs_per_line; ++d)
    *next++ = dof_index(d);
}



template <>
inline
unsigned int
DoFObjectAccessor<1,DoFHandler<2> >::dof_index (const unsigned int i,
						const unsigned int fe_index) const
{
  Assert (fe_index == DoFHandler<2>::default_fe_index,
	  ExcMessage ("Only the default FE index is allowed for non-hp DoFHandler objects"));
  Assert (static_cast<unsigned int>(this->present_level) < this->dof_handler->levels.size(),
          ExcMessage ("DoFHandler not initialized"));

  Assert (this->dof_handler != 0, DoFAccessor<DoFHandler<2> >::ExcInvalidObject());
				   // make sure a FE has been selected
				   // and enough room was reserved
  Assert (&this->dof_handler->get_fe() != 0, DoFAccessor<DoFHandler<2> >::ExcInvalidObject());
  Assert (i<this->dof_handler->get_fe().dofs_per_line,
	  ExcIndexRange (i, 0, this->dof_handler->get_fe().dofs_per_line));

  return this->dof_handler->levels[this->present_level]
    ->line_dofs[this->present_index*this->dof_handler->get_fe().dofs_per_line+i];
}



template <>
inline
unsigned int
DoFObjectAccessor<1,DoFHandler<2> >::vertex_dof_index (const unsigned int vertex,
						       const unsigned int i,
						       const unsigned int fe_index) const
{
  Assert (fe_index == DoFHandler<2>::default_fe_index,
	  ExcMessage ("Only the default FE index is allowed for non-hp DoFHandler objects"));
  Assert (static_cast<unsigned int>(this->present_level) < this->dof_handler->levels.size(),
          ExcMessage ("DoFHandler not initialized"));

  Assert (this->dof_handler != 0, DoFAccessor<DoFHandler<2> >::ExcInvalidObject());
  Assert (&this->dof_handler->get_fe() != 0, DoFAccessor<DoFHandler<2> >::ExcInvalidObject());
  Assert (vertex<2, ExcIndexRange (i,0,2));
  Assert (i<this->dof_handler->get_fe().dofs_per_vertex,
	  ExcIndexRange (i, 0, this->dof_handler->get_fe().dofs_per_vertex));

  const unsigned int dof_number = (this->vertex_index(vertex) *
				   this->dof_handler->get_fe().dofs_per_vertex +
				   i);
  return this->dof_handler->vertex_dofs[dof_number];
}



template <>
inline
void
DoFObjectAccessor<1,DoFHandler<2> >::get_dof_indices (std::vector<unsigned int> &dof_indices) const
{
  Assert (static_cast<unsigned int>(this->present_level) < this->dof_handler->levels.size(),
          ExcMessage ("DoFHandler not initialized"));

  Assert (this->dof_handler != 0, DoFAccessor<DoFHandler<2> >::ExcInvalidObject());
  Assert (&this->dof_handler->get_fe() != 0, DoFAccessor<DoFHandler<2> >::ExcInvalidObject());
  Assert (dof_indices.size() == (2*this->dof_handler->get_fe().dofs_per_vertex +
				 this->dof_handler->get_fe().dofs_per_line),
	  DoFAccessor<DoFHandler<2> >::ExcVectorDoesNotMatch());

				   // this function really only makes
				   // sense on non-active objects if
				   // all degrees of freedom are
				   // located on vertices, since
				   // otherwise there are degrees of
				   // freedom on sub-objects which are
				   // not allocated for this
				   // non-active thing
  Assert (!this->has_children() ||
	  (this->dof_handler->get_fe().dofs_per_cell ==
	   2*this->dof_handler->get_fe().dofs_per_vertex),
	  DoFAccessor<DoFHandler<2> >::ExcNotActive());
  
  const unsigned int dofs_per_vertex = this->dof_handler->get_fe().dofs_per_vertex,
		     dofs_per_line   = this->dof_handler->get_fe().dofs_per_line;
  std::vector<unsigned int>::iterator next = dof_indices.begin();
  for (unsigned int vertex=0; vertex<2; ++vertex)
    for (unsigned int d=0; d<dofs_per_vertex; ++d)
      *next++ = vertex_dof_index(vertex,d);
  for (unsigned int d=0; d<dofs_per_line; ++d)
    *next++ = dof_index(d);
}



template <>
inline
unsigned int
DoFObjectAccessor<1,DoFHandler<3> >::dof_index (const unsigned int i,
						const unsigned int fe_index) const
{
  Assert (fe_index == DoFHandler<3>::default_fe_index,
	  ExcMessage ("Only the default FE index is allowed for non-hp DoFHandler objects"));
  Assert (static_cast<unsigned int>(this->present_level) < this->dof_handler->levels.size(),
          ExcMessage ("DoFHandler not initialized"));

  Assert (this->dof_handler != 0, DoFAccessor<DoFHandler<3> >::ExcInvalidObject());
				   // make sure a FE has been selected
				   // and enough room was reserved
  Assert (&this->dof_handler->get_fe() != 0, DoFAccessor<DoFHandler<3> >::ExcInvalidObject());
  Assert (i<this->dof_handler->get_fe().dofs_per_line,
	  ExcIndexRange (i, 0, this->dof_handler->get_fe().dofs_per_line));

  return this->dof_handler->levels[this->present_level]
    ->line_dofs[this->present_index*this->dof_handler->get_fe().dofs_per_line+i];
}



template <>
inline
unsigned int
DoFObjectAccessor<1,DoFHandler<3> >::vertex_dof_index (const unsigned int vertex,
						       const unsigned int i,
						       const unsigned int fe_index) const
{
  Assert (fe_index == DoFHandler<3>::default_fe_index,
	  ExcMessage ("Only the default FE index is allowed for non-hp DoFHandler objects"));
  Assert (static_cast<unsigned int>(this->present_level) < this->dof_handler->levels.size(),
          ExcMessage ("DoFHandler not initialized"));

  Assert (this->dof_handler != 0, DoFAccessor<DoFHandler<3> >::ExcInvalidObject());
  Assert (&this->dof_handler->get_fe() != 0, DoFAccessor<DoFHandler<3> >::ExcInvalidObject());
  Assert (vertex<2, ExcIndexRange (i,0,2));
  Assert (i<this->dof_handler->get_fe().dofs_per_vertex,
	  ExcIndexRange (i, 0, this->dof_handler->get_fe().dofs_per_vertex));

  const unsigned int dof_number = (this->vertex_index(vertex) *
				   this->dof_handler->get_fe().dofs_per_vertex +
				   i);
  return this->dof_handler->vertex_dofs[dof_number];
}



template <>
inline
void
DoFObjectAccessor<1,DoFHandler<3> >::get_dof_indices (std::vector<unsigned int> &dof_indices) const
{
  Assert (static_cast<unsigned int>(this->present_level) < this->dof_handler->levels.size(),
          ExcMessage ("DoFHandler not initialized"));

  Assert (this->dof_handler != 0, DoFAccessor<DoFHandler<3> >::ExcInvalidObject());
  Assert (&this->dof_handler->get_fe() != 0, DoFAccessor<DoFHandler<3> >::ExcInvalidObject());
  Assert (dof_indices.size() == (2*this->dof_handler->get_fe().dofs_per_vertex +
				 this->dof_handler->get_fe().dofs_per_line),
	  DoFAccessor<DoFHandler<3> >::ExcVectorDoesNotMatch());

				   // this function really only makes
				   // sense on non-active objects if
				   // all degrees of freedom are
				   // located on vertices, since
				   // otherwise there are degrees of
				   // freedom on sub-objects which are
				   // not allocated for this
				   // non-active thing
  Assert (!this->has_children() ||
	  (this->dof_handler->get_fe().dofs_per_cell ==
	   2*this->dof_handler->get_fe().dofs_per_vertex),
	  DoFAccessor<DoFHandler<3> >::ExcNotActive());
  
  const unsigned int dofs_per_vertex = this->dof_handler->get_fe().dofs_per_vertex,
		     dofs_per_line   = this->dof_handler->get_fe().dofs_per_line;
  std::vector<unsigned int>::iterator next = dof_indices.begin();
  for (unsigned int vertex=0; vertex<2; ++vertex)
    for (unsigned int d=0; d<dofs_per_vertex; ++d)
      *next++ = vertex_dof_index(vertex,d);
  for (unsigned int d=0; d<dofs_per_line; ++d)
    *next++ = dof_index(d);
}



/*--------------- Functions: DoFObjectAccessor<2,DoFHandler> -----------*/

template <>
inline
unsigned int
DoFObjectAccessor<2,DoFHandler<2> >::dof_index (const unsigned int i,
						const unsigned int fe_index) const
{
  Assert (fe_index == DoFHandler<2>::default_fe_index,
	  ExcMessage ("Only the default FE index is allowed for non-hp DoFHandler objects"));
  Assert (this->dof_handler != 0,
	  DoFAccessor<DoFHandler<2> >::ExcInvalidObject());
				   // make sure a FE has been selected
				   // and enough room was reserved
  Assert (&this->dof_handler->get_fe() != 0,
	  DoFAccessor<DoFHandler<2> >::ExcInvalidObject());
  Assert (i<this->dof_handler->get_fe().dofs_per_quad,
	  ExcIndexRange (i, 0, this->dof_handler->get_fe().dofs_per_quad));
  Assert (static_cast<unsigned int>(this->present_level) < this->dof_handler->levels.size(),
          ExcMessage ("DoFHandler not initialized"));

  return this->dof_handler->levels[this->present_level]
    ->quad_dofs[this->present_index*this->dof_handler->get_fe().dofs_per_quad+i];
}



template <>
inline
unsigned int
DoFObjectAccessor<2,DoFHandler<2> >::vertex_dof_index (const unsigned int vertex,
						       const unsigned int i,
						       const unsigned int fe_index) const
{
  Assert (fe_index == DoFHandler<2>::default_fe_index,
	  ExcMessage ("Only the default FE index is allowed for non-hp DoFHandler objects"));
  Assert (this->dof_handler != 0,
	  DoFAccessor<DoFHandler<2> >::ExcInvalidObject());
  Assert (&this->dof_handler->get_fe() != 0,
	  DoFAccessor<DoFHandler<2> >::ExcInvalidObject());
  Assert (vertex<4, ExcIndexRange (i,0,4));
  Assert (i<this->dof_handler->get_fe().dofs_per_vertex,
	  ExcIndexRange (i, 0, this->dof_handler->get_fe().dofs_per_vertex));
  Assert (static_cast<unsigned int>(this->present_level) < this->dof_handler->levels.size(),
          ExcMessage ("DoFHandler not initialized"));

  const unsigned int dof_number = (this->vertex_index(vertex) *
				   this->dof_handler->get_fe().dofs_per_vertex +
				   i);
  return this->dof_handler->vertex_dofs[dof_number];
}



template <>
inline
void
DoFObjectAccessor<2,DoFHandler<2> >::get_dof_indices (std::vector<unsigned int> &dof_indices) const
{
  Assert (this->dof_handler != 0,
	  DoFAccessor<DoFHandler<2> >::ExcInvalidObject());
  Assert (&this->dof_handler->get_fe() != 0,
	  DoFAccessor<DoFHandler<2> >::ExcInvalidObject());
  Assert (dof_indices.size() == (4*this->dof_handler->get_fe().dofs_per_vertex +
				 4*this->dof_handler->get_fe().dofs_per_line +
				 this->dof_handler->get_fe().dofs_per_quad),
	  DoFAccessor<DoFHandler<2> >::ExcVectorDoesNotMatch());

				   // this function really only makes
				   // sense on non-active objects if
				   // all degrees of freedom are
				   // located on vertices, since
				   // otherwise there are degrees of
				   // freedom on sub-objects which are
				   // not allocated for this
				   // non-active thing
  Assert (!this->has_children() ||
	  (this->dof_handler->get_fe().dofs_per_cell ==
	   4*this->dof_handler->get_fe().dofs_per_vertex),
	  DoFAccessor<DoFHandler<2> >::ExcNotActive());
  
  Assert (static_cast<unsigned int>(this->present_level) < this->dof_handler->levels.size(),
          ExcMessage ("DoFHandler not initialized"));
  
  const unsigned int dofs_per_vertex = this->dof_handler->get_fe().dofs_per_vertex,
		     dofs_per_line   = this->dof_handler->get_fe().dofs_per_line,
		     dofs_per_quad   = this->dof_handler->get_fe().dofs_per_quad;
  std::vector<unsigned int>::iterator next = dof_indices.begin();
  for (unsigned int vertex=0; vertex<4; ++vertex)
    for (unsigned int d=0; d<dofs_per_vertex; ++d)
      *next++ = vertex_dof_index(vertex,d);
  for (unsigned int line=0; line<4; ++line)
    for (unsigned int d=0; d<dofs_per_line; ++d)
      *next++ = this->line(line)->dof_index(d);
  for (unsigned int d=0; d<dofs_per_quad; ++d)
    *next++ = dof_index(d);
}



template <>
inline
unsigned int
DoFObjectAccessor<2,DoFHandler<3> >::dof_index (const unsigned int i,
						const unsigned int fe_index) const
{
  Assert (fe_index == DoFHandler<3>::default_fe_index,
	  ExcMessage ("Only the default FE index is allowed for non-hp DoFHandler objects"));
  Assert (this->dof_handler != 0,
	  DoFAccessor<DoFHandler<3> >::ExcInvalidObject());
				   // make sure a FE has been selected
				   // and enough room was reserved
  Assert (&this->dof_handler->get_fe() != 0,
	  DoFAccessor<DoFHandler<3> >::ExcInvalidObject());
  Assert (i<this->dof_handler->get_fe().dofs_per_quad,
	  ExcIndexRange (i, 0, this->dof_handler->get_fe().dofs_per_quad));
  Assert (static_cast<unsigned int>(this->present_level) < this->dof_handler->levels.size(),
          ExcMessage ("DoFHandler not initialized"));

  return this->dof_handler->levels[this->present_level]
    ->quad_dofs[this->present_index*this->dof_handler->get_fe().dofs_per_quad+i];
}



template <>
inline
unsigned int
DoFObjectAccessor<2,DoFHandler<3> >::vertex_dof_index (const unsigned int vertex,
						       const unsigned int i,
						       const unsigned int fe_index) const
{
  Assert (fe_index == DoFHandler<3>::default_fe_index,
	  ExcMessage ("Only the default FE index is allowed for non-hp DoFHandler objects"));
  Assert (this->dof_handler != 0,
	  DoFAccessor<DoFHandler<3> >::ExcInvalidObject());
  Assert (&this->dof_handler->get_fe() != 0,
	  DoFAccessor<DoFHandler<3> >::ExcInvalidObject());
  Assert (vertex<4, ExcIndexRange (i,0,4));
  Assert (i<this->dof_handler->get_fe().dofs_per_vertex,
	  ExcIndexRange (i, 0, this->dof_handler->get_fe().dofs_per_vertex));
  Assert (static_cast<unsigned int>(this->present_level) < this->dof_handler->levels.size(),
          ExcMessage ("DoFHandler not initialized"));

  const unsigned int dof_number = (this->vertex_index(vertex) *
				   this->dof_handler->get_fe().dofs_per_vertex +
				   i);
  return this->dof_handler->vertex_dofs[dof_number];
}



template <>
inline
void
DoFObjectAccessor<2,DoFHandler<3> >::get_dof_indices (std::vector<unsigned int> &dof_indices) const
{
  Assert (this->dof_handler != 0,
	  DoFAccessor<DoFHandler<3> >::ExcInvalidObject());
  Assert (&this->dof_handler->get_fe() != 0,
	  DoFAccessor<DoFHandler<3> >::ExcInvalidObject());
  Assert (dof_indices.size() == (4*this->dof_handler->get_fe().dofs_per_vertex +
				 4*this->dof_handler->get_fe().dofs_per_line +
				 this->dof_handler->get_fe().dofs_per_quad),
	  DoFAccessor<DoFHandler<3> >::ExcVectorDoesNotMatch());

				   // this function really only makes
				   // sense on non-active objects if
				   // all degrees of freedom are
				   // located on vertices, since
				   // otherwise there are degrees of
				   // freedom on sub-objects which are
				   // not allocated for this
				   // non-active thing
  Assert (!this->has_children() ||
	  (this->dof_handler->get_fe().dofs_per_cell ==
	   4*this->dof_handler->get_fe().dofs_per_vertex),
	  DoFAccessor<DoFHandler<3> >::ExcNotActive());
  
  Assert (static_cast<unsigned int>(this->present_level) < this->dof_handler->levels.size(),
          ExcMessage ("DoFHandler not initialized"));
  
  const unsigned int dofs_per_vertex = this->dof_handler->get_fe().dofs_per_vertex,
		     dofs_per_line   = this->dof_handler->get_fe().dofs_per_line,
		     dofs_per_quad   = this->dof_handler->get_fe().dofs_per_quad;
  std::vector<unsigned int>::iterator next = dof_indices.begin();
  for (unsigned int vertex=0; vertex<4; ++vertex)
    for (unsigned int d=0; d<dofs_per_vertex; ++d)
      *next++ = vertex_dof_index(vertex,d);
  for (unsigned int line=0; line<4; ++line)
    for (unsigned int d=0; d<dofs_per_line; ++d)
      *next++ = this->line(line)->dof_index(d);
  for (unsigned int d=0; d<dofs_per_quad; ++d)
    *next++ = dof_index(d);
}



/*--------------- Functions: DoFObjectAccessor<3,DoFHandler> -----------*/




template <>
inline
unsigned int
DoFObjectAccessor<3,DoFHandler<3> >::dof_index (const unsigned int i,
						const unsigned int fe_index) const
{
  Assert (fe_index == DoFHandler<3>::default_fe_index,
	  ExcMessage ("Only the default FE index is allowed for non-hp DoFHandler objects"));
  Assert (this->dof_handler != 0,
	  DoFAccessor<DoFHandler<3> >::ExcInvalidObject());
				   // make sure a FE has been selected
				   // and enough room was reserved
  Assert (&this->dof_handler->get_fe() != 0,
	  DoFAccessor<DoFHandler<3> >::ExcInvalidObject());
  Assert (i<this->dof_handler->get_fe().dofs_per_hex,
	  ExcIndexRange (i, 0, this->dof_handler->get_fe().dofs_per_hex));
  Assert (static_cast<unsigned int>(this->present_level) < this->dof_handler->levels.size(),
          ExcMessage ("DoFHandler not initialized"));

  return this->dof_handler->levels[this->present_level]
    ->hex_dofs[this->present_index*this->dof_handler->get_fe().dofs_per_hex+i];
}



template <>
inline
unsigned int
DoFObjectAccessor<3,DoFHandler<3> >::vertex_dof_index (const unsigned int vertex,
						       const unsigned int i,
						       const unsigned int fe_index) const
{
  Assert (fe_index == DoFHandler<3>::default_fe_index,
	  ExcMessage ("Only the default FE index is allowed for non-hp DoFHandler objects"));
  Assert (this->dof_handler != 0,
	  DoFAccessor<DoFHandler<3> >::ExcInvalidObject());
  Assert (&this->dof_handler->get_fe() != 0,
	  DoFAccessor<DoFHandler<3> >::ExcInvalidObject());
  Assert (vertex<8, ExcIndexRange (i,0,8));
  Assert (i<this->dof_handler->get_fe().dofs_per_vertex,
	  ExcIndexRange (i, 0, this->dof_handler->get_fe().dofs_per_vertex));
  Assert (static_cast<unsigned int>(this->present_level) < this->dof_handler->levels.size(),
          ExcMessage ("DoFHandler not initialized"));

  const unsigned int dof_number = (this->vertex_index(vertex) *
				   this->dof_handler->get_fe().dofs_per_vertex +
				   i);
  return this->dof_handler->vertex_dofs[dof_number];
}



template <>
inline
void
DoFObjectAccessor<3,DoFHandler<3> >::get_dof_indices (std::vector<unsigned int> &dof_indices) const
{
  Assert (this->dof_handler != 0,
	  DoFAccessor<DoFHandler<3> >::ExcInvalidObject());
  Assert (&this->dof_handler->get_fe() != 0,
	  DoFAccessor<DoFHandler<3> >::ExcInvalidObject());
  Assert (dof_indices.size() == (8*this->dof_handler->get_fe().dofs_per_vertex +
				 12*this->dof_handler->get_fe().dofs_per_line +
				 6*this->dof_handler->get_fe().dofs_per_quad +
				 this->dof_handler->get_fe().dofs_per_hex),
	  DoFAccessor<DoFHandler<3> >::ExcVectorDoesNotMatch());
  Assert (static_cast<unsigned int>(this->present_level) < this->dof_handler->levels.size(),
          ExcMessage ("DoFHandler not initialized"));

				   // this function really only makes
				   // sense on non-active objects if
				   // all degrees of freedom are
				   // located on vertices, since
				   // otherwise there are degrees of
				   // freedom on sub-objects which are
				   // not allocated for this
				   // non-active thing
  Assert (!this->has_children() ||
	  (this->dof_handler->get_fe().dofs_per_cell ==
	   8*this->dof_handler->get_fe().dofs_per_vertex),
	  DoFAccessor<DoFHandler<3> >::ExcNotActive());
  Assert (static_cast<unsigned int>(this->present_level) < this->dof_handler->levels.size(),
          ExcMessage ("DoFHandler not initialized"));
  
  const unsigned int dofs_per_vertex = this->dof_handler->get_fe().dofs_per_vertex,
		     dofs_per_line   = this->dof_handler->get_fe().dofs_per_line,
		     dofs_per_quad   = this->dof_handler->get_fe().dofs_per_quad,
		     dofs_per_hex    = this->dof_handler->get_fe().dofs_per_hex;
  std::vector<unsigned int>::iterator next = dof_indices.begin();
  for (unsigned int vertex=0; vertex<8; ++vertex)
    for (unsigned int d=0; d<dofs_per_vertex; ++d)
      *next++ = vertex_dof_index(vertex,d);
  for (unsigned int line=0; line<12; ++line)
    for (unsigned int d=0; d<dofs_per_line; ++d)
      *next++ = this->line(line)->dof_index(d);
  for (unsigned int quad=0; quad<6; ++quad)
    for (unsigned int d=0; d<dofs_per_quad; ++d)
      *next++ = this->quad(quad)->dof_index(d);
  for (unsigned int d=0; d<dofs_per_hex; ++d)
    *next++ = dof_index(d);
}



/* -------------- hp vertex dofs stuff ------------------------------------ */

template <>
inline
unsigned int
DoFObjectAccessor<1,hp::DoFHandler<1> >::vertex_dof_index (const unsigned int vertex,
							   const unsigned int i,
							   const unsigned int fe_index) const
{
  return internal::hp::DoFLevel<0>::
    get_vertex_dof_index (*this->dof_handler,
                          this->vertex_index(vertex),
                          fe_index,
                          i);
}



template <>
inline
unsigned int
DoFObjectAccessor<1,hp::DoFHandler<2> >::vertex_dof_index (const unsigned int vertex,
							   const unsigned int i,
							   const unsigned int fe_index) const
{
  return internal::hp::DoFLevel<0>::
    get_vertex_dof_index (*this->dof_handler,
                          this->vertex_index(vertex),
                          fe_index,
                          i);
}


template <>
inline
unsigned int
DoFObjectAccessor<1,hp::DoFHandler<3> >::vertex_dof_index (const unsigned int vertex,
							   const unsigned int i,
							   const unsigned int fe_index) const
{
  return internal::hp::DoFLevel<0>::
    get_vertex_dof_index (*this->dof_handler,
                          this->vertex_index(vertex),
                          fe_index,
                          i);
}


template <>
inline
unsigned int
DoFObjectAccessor<2,hp::DoFHandler<2> >::vertex_dof_index (const unsigned int vertex,
							   const unsigned int i,
							   const unsigned int fe_index) const
{
  return internal::hp::DoFLevel<0>::
    get_vertex_dof_index (*this->dof_handler,
                          this->vertex_index(vertex),
                          fe_index,
                          i);
}



template <>
inline
unsigned int
DoFObjectAccessor<2,hp::DoFHandler<3> >::vertex_dof_index (const unsigned int vertex,
							   const unsigned int i,
							   const unsigned int fe_index) const
{
  return internal::hp::DoFLevel<0>::
    get_vertex_dof_index (*this->dof_handler,
                          this->vertex_index(vertex),
                          fe_index,
                          i);
}



template <>
inline
unsigned int
DoFObjectAccessor<3,hp::DoFHandler<3> >::vertex_dof_index (const unsigned int vertex,
							   const unsigned int i,
							   const unsigned int fe_index) const
{
  return internal::hp::DoFLevel<0>::
    get_vertex_dof_index (*this->dof_handler,
                          this->vertex_index(vertex),
                          fe_index,
                          i);
}


/*--------------- Functions: DoFObjectAccessor<1,hp::DoFHandler> -----------*/

template <>
inline
unsigned int
DoFObjectAccessor<1,hp::DoFHandler<1> >::dof_index (const unsigned int i,
						    const unsigned int fe_index) const
{
  return this->dof_handler->levels[this->present_level]
    ->get_line_dof_index (*this->dof_handler,
                          this->present_index,
                          fe_index,
                          i);
}



template <>
inline
void
DoFObjectAccessor<1,hp::DoFHandler<1> >::get_dof_indices (std::vector<unsigned int> &dof_indices) const
{
  Assert (static_cast<unsigned int>(this->present_level) < this->dof_handler->levels.size(),
          ExcMessage ("DoFHandler not initialized"));

  Assert (this->dof_handler != 0, DoFAccessor<hp::DoFHandler<1> >::ExcInvalidObject());
  Assert (&this->dof_handler->get_fe() != 0, DoFAccessor<hp::DoFHandler<1> >::ExcInvalidObject());
  Assert (dof_indices.size() == (2*this->dof_handler->get_fe()[0].dofs_per_vertex +
				 this->dof_handler->get_fe()[0].dofs_per_line),
	  DoFAccessor<hp::DoFHandler<1> >::ExcVectorDoesNotMatch());

				   // this function really only makes
				   // sense on non-active objects if
				   // all degrees of freedom are
				   // located on vertices, since
				   // otherwise there are degrees of
				   // freedom on sub-objects which are
				   // not allocated for this
				   // non-active thing
  Assert (!this->has_children() ||
	  (this->dof_handler->get_fe()[0].dofs_per_cell ==
	   2*this->dof_handler->get_fe()[0].dofs_per_vertex),
	  DoFAccessor<hp::DoFHandler<1> >::ExcNotActive());
  
  const unsigned int dofs_per_vertex = this->dof_handler->get_fe()[0].dofs_per_vertex,
		     dofs_per_line   = this->dof_handler->get_fe()[0].dofs_per_line;
  std::vector<unsigned int>::iterator next = dof_indices.begin();
  for (unsigned int vertex=0; vertex<2; ++vertex)
    for (unsigned int d=0; d<dofs_per_vertex; ++d)
      *next++ = vertex_dof_index(vertex,d);
  for (unsigned int d=0; d<dofs_per_line; ++d)
    *next++ = dof_index(d);
}



template <>
inline
unsigned int
DoFObjectAccessor<1,hp::DoFHandler<2> >::dof_index (const unsigned int i,
						    const unsigned int fe_index) const
{
  return this->dof_handler->levels[this->present_level]
    ->get_line_dof_index (*this->dof_handler,
                          this->present_index,
                          fe_index,
                          i);
}



template <>
inline
void
DoFObjectAccessor<1,hp::DoFHandler<2> >::get_dof_indices (std::vector<unsigned int> &dof_indices) const
{
  Assert (static_cast<unsigned int>(this->present_level) < this->dof_handler->levels.size(),
          ExcMessage ("DoFHandler not initialized"));

  Assert (this->dof_handler != 0, DoFAccessor<hp::DoFHandler<2> >::ExcInvalidObject());
  Assert (&this->dof_handler->get_fe() != 0, DoFAccessor<hp::DoFHandler<2> >::ExcInvalidObject());
  Assert (dof_indices.size() == (2*this->dof_handler->get_fe()[0].dofs_per_vertex +
				 this->dof_handler->get_fe()[0].dofs_per_line),
	  DoFAccessor<hp::DoFHandler<2> >::ExcVectorDoesNotMatch());

				   // this function really only makes
				   // sense on non-active objects if
				   // all degrees of freedom are
				   // located on vertices, since
				   // otherwise there are degrees of
				   // freedom on sub-objects which are
				   // not allocated for this
				   // non-active thing
  Assert (!this->has_children() ||
	  (this->dof_handler->get_fe()[0].dofs_per_cell ==
	   2*this->dof_handler->get_fe()[0].dofs_per_vertex),
	  DoFAccessor<hp::DoFHandler<2> >::ExcNotActive());
  
  const unsigned int dofs_per_vertex = this->dof_handler->get_fe()[0].dofs_per_vertex,
		     dofs_per_line   = this->dof_handler->get_fe()[0].dofs_per_line;
  std::vector<unsigned int>::iterator next = dof_indices.begin();
  for (unsigned int vertex=0; vertex<2; ++vertex)
    for (unsigned int d=0; d<dofs_per_vertex; ++d)
      *next++ = vertex_dof_index(vertex,d);
  for (unsigned int d=0; d<dofs_per_line; ++d)
    *next++ = dof_index(d);
}



template <>
inline
unsigned int
DoFObjectAccessor<1,hp::DoFHandler<3> >::dof_index (const unsigned int i,
						    const unsigned int fe_index) const
{
  return this->dof_handler->levels[this->present_level]
    ->get_line_dof_index (*this->dof_handler,
                          this->present_index,
                          fe_index,
                          i);
}



template <>
inline
void
DoFObjectAccessor<1,hp::DoFHandler<3> >::get_dof_indices (std::vector<unsigned int> &dof_indices) const
{
  Assert (static_cast<unsigned int>(this->present_level) < this->dof_handler->levels.size(),
          ExcMessage ("DoFHandler not initialized"));

  Assert (this->dof_handler != 0, DoFAccessor<hp::DoFHandler<3> >::ExcInvalidObject());
  Assert (&this->dof_handler->get_fe() != 0, DoFAccessor<hp::DoFHandler<3> >::ExcInvalidObject());
  Assert (dof_indices.size() == (2*this->dof_handler->get_fe()[0].dofs_per_vertex +
				 this->dof_handler->get_fe()[0].dofs_per_line),
	  DoFAccessor<hp::DoFHandler<3> >::ExcVectorDoesNotMatch());

				   // this function really only makes
				   // sense on non-active objects if
				   // all degrees of freedom are
				   // located on vertices, since
				   // otherwise there are degrees of
				   // freedom on sub-objects which are
				   // not allocated for this
				   // non-active thing
  Assert (!this->has_children() ||
	  (this->dof_handler->get_fe()[0].dofs_per_cell ==
	   2*this->dof_handler->get_fe()[0].dofs_per_vertex),
	  DoFAccessor<hp::DoFHandler<3> >::ExcNotActive());
  
  const unsigned int dofs_per_vertex = this->dof_handler->get_fe()[0].dofs_per_vertex,
		     dofs_per_line   = this->dof_handler->get_fe()[0].dofs_per_line;
  std::vector<unsigned int>::iterator next = dof_indices.begin();
  for (unsigned int vertex=0; vertex<2; ++vertex)
    for (unsigned int d=0; d<dofs_per_vertex; ++d)
      *next++ = vertex_dof_index(vertex,d);
  for (unsigned int d=0; d<dofs_per_line; ++d)
    *next++ = dof_index(d);
}



template <>
inline
void
DoFObjectAccessor<1,hp::DoFHandler<1> >::set_dof_index (const unsigned int i,
							const unsigned int index,
							const unsigned int fe_index) const
{
  this->dof_handler->levels[this->present_level]
    ->set_line_dof_index (*this->dof_handler,
                          this->present_index,
                          fe_index,
                          i,
                          index);
}



template <>
inline
void
DoFObjectAccessor<1, hp::DoFHandler<2> >::set_dof_index (const unsigned int i,
							 const unsigned int index,
							 const unsigned int fe_index) const
{
  this->dof_handler->levels[this->present_level]
    ->set_line_dof_index (*this->dof_handler,
                          this->present_index,
                          fe_index,
                          i,
                          index);
}



template <>
inline
void
DoFObjectAccessor<1, hp::DoFHandler<3> >::set_dof_index (const unsigned int i,
							 const unsigned int index,
							 const unsigned int fe_index) const
{
  this->dof_handler->levels[this->present_level]
    ->set_line_dof_index (*this->dof_handler,
                          this->present_index,
                          fe_index,
                          i,
                          index);
}




/*------------- Functions: DoFObjectAccessor<2,hp::DoFHandler> -------------*/

template <>
inline
unsigned int DoFObjectAccessor<2,hp::DoFHandler<2> >::dof_index (const unsigned int i,
								 const unsigned int fe_index) const
{
  return this->dof_handler->levels[this->present_level]
    ->get_quad_dof_index (*this->dof_handler,
                          this->present_index,
                          fe_index,
                          i);
}



template <>
inline
unsigned int DoFObjectAccessor<2,hp::DoFHandler<3> >::dof_index (const unsigned int i,
								 const unsigned int fe_index) const
{
  return this->dof_handler->levels[this->present_level]
    ->get_quad_dof_index (*this->dof_handler,
                          this->present_index,
                          fe_index,
                          i);
}



template <>
inline
void
DoFObjectAccessor<2, hp::DoFHandler<2> >::set_dof_index (const unsigned int i,
							 const unsigned int index,
							 const unsigned int fe_index) const
{
  this->dof_handler->levels[this->present_level]
    ->set_quad_dof_index (*this->dof_handler,
                          this->present_index,
                          fe_index,
                          i,
                          index);
}



template <>
inline
void
DoFObjectAccessor<2, hp::DoFHandler<3> >::set_dof_index (const unsigned int i,
							 const unsigned int index,
							 const unsigned int fe_index) const
{
  this->dof_handler->levels[this->present_level]
    ->set_quad_dof_index (*this->dof_handler,
                          this->present_index,
                          fe_index,
                          i,
                          index);
}



template <>
inline
void
DoFObjectAccessor<2,hp::DoFHandler<2> >::get_dof_indices (std::vector<unsigned int> &dof_indices) const
{
  Assert (this->dof_handler != 0,
	  DoFAccessor<hp::DoFHandler<2> >::ExcInvalidObject());
  Assert (&this->dof_handler->get_fe() != 0,
	  DoFAccessor<hp::DoFHandler<2> >::ExcInvalidObject());
  Assert (dof_indices.size() == (4*this->dof_handler->get_fe()[0].dofs_per_vertex +
				 4*this->dof_handler->get_fe()[0].dofs_per_line +
				 this->dof_handler->get_fe()[0].dofs_per_quad),
	  DoFAccessor<hp::DoFHandler<2> >::ExcVectorDoesNotMatch());

				   // this function really only makes
				   // sense on non-active objects if
				   // all degrees of freedom are
				   // located on vertices, since
				   // otherwise there are degrees of
				   // freedom on sub-objects which are
				   // not allocated for this
				   // non-active thing
  Assert (!this->has_children() ||
	  (this->dof_handler->get_fe()[0].dofs_per_cell ==
	   4*this->dof_handler->get_fe()[0].dofs_per_vertex),
	  DoFAccessor<hp::DoFHandler<2> >::ExcNotActive());
  
  Assert (static_cast<unsigned int>(this->present_level) < this->dof_handler->levels.size(),
          ExcMessage ("DoFHandler not initialized"));
  
  const unsigned int dofs_per_vertex = this->dof_handler->get_fe()[0].dofs_per_vertex,
		     dofs_per_line   = this->dof_handler->get_fe()[0].dofs_per_line,
		     dofs_per_quad   = this->dof_handler->get_fe()[0].dofs_per_quad;
  std::vector<unsigned int>::iterator next = dof_indices.begin();
  for (unsigned int vertex=0; vertex<4; ++vertex)
    for (unsigned int d=0; d<dofs_per_vertex; ++d)
      *next++ = vertex_dof_index(vertex,d);
  for (unsigned int line=0; line<4; ++line)
    for (unsigned int d=0; d<dofs_per_line; ++d)
      *next++ = this->line(line)->dof_index(d);
  for (unsigned int d=0; d<dofs_per_quad; ++d)
    *next++ = dof_index(d);
}



template <>
inline
void
DoFObjectAccessor<2,hp::DoFHandler<3> >::get_dof_indices (std::vector<unsigned int> &dof_indices) const
{
  Assert (this->dof_handler != 0,
	  DoFAccessor<hp::DoFHandler<3> >::ExcInvalidObject());
  Assert (&this->dof_handler->get_fe() != 0,
	  DoFAccessor<hp::DoFHandler<3> >::ExcInvalidObject());
  Assert (dof_indices.size() == (4*this->dof_handler->get_fe()[0].dofs_per_vertex +
				 4*this->dof_handler->get_fe()[0].dofs_per_line +
				 this->dof_handler->get_fe()[0].dofs_per_quad),
	  DoFAccessor<hp::DoFHandler<3> >::ExcVectorDoesNotMatch());

				   // this function really only makes
				   // sense on non-active objects if
				   // all degrees of freedom are
				   // located on vertices, since
				   // otherwise there are degrees of
				   // freedom on sub-objects which are
				   // not allocated for this
				   // non-active thing
  Assert (!this->has_children() ||
	  (this->dof_handler->get_fe()[0].dofs_per_cell ==
	   4*this->dof_handler->get_fe()[0].dofs_per_vertex),
	  DoFAccessor<hp::DoFHandler<3> >::ExcNotActive());
  
  Assert (static_cast<unsigned int>(this->present_level) < this->dof_handler->levels.size(),
          ExcMessage ("DoFHandler not initialized"));
  
  const unsigned int dofs_per_vertex = this->dof_handler->get_fe()[0].dofs_per_vertex,
		     dofs_per_line   = this->dof_handler->get_fe()[0].dofs_per_line,
		     dofs_per_quad   = this->dof_handler->get_fe()[0].dofs_per_quad;
  std::vector<unsigned int>::iterator next = dof_indices.begin();
  for (unsigned int vertex=0; vertex<4; ++vertex)
    for (unsigned int d=0; d<dofs_per_vertex; ++d)
      *next++ = vertex_dof_index(vertex,d);
  for (unsigned int line=0; line<4; ++line)
    for (unsigned int d=0; d<dofs_per_line; ++d)
      *next++ = this->line(line)->dof_index(d);
  for (unsigned int d=0; d<dofs_per_quad; ++d)
    *next++ = dof_index(d);
}





/*------------- Functions: DoFObjectAccessor<3,hp::DoFHandler> -------------*/

template <>
inline
unsigned int
DoFObjectAccessor<3,hp::DoFHandler<3> >::dof_index (const unsigned int i,
						    const unsigned int fe_index) const
{
  return this->dof_handler->levels[this->present_level]
    ->get_hex_dof_index (*this->dof_handler,
                         this->present_index,
                         fe_index,
                         i);
}



template <>
inline
void
DoFObjectAccessor<3, hp::DoFHandler<3> >::set_dof_index (const unsigned int i,
							 const unsigned int index,
							 const unsigned int fe_index) const
{
  this->dof_handler->levels[this->present_level]
    ->set_hex_dof_index (*this->dof_handler,
                         this->present_index,
                         fe_index,
                         i,
                         index);
}



template <>
inline
void
DoFObjectAccessor<3,hp::DoFHandler<3> >::get_dof_indices (std::vector<unsigned int> &dof_indices) const
{
  Assert (this->dof_handler != 0,
	  DoFAccessor<hp::DoFHandler<3> >::ExcInvalidObject());
  Assert (&this->dof_handler->get_fe() != 0,
	  DoFAccessor<hp::DoFHandler<3> >::ExcInvalidObject());
  Assert (dof_indices.size() == (8*this->dof_handler->get_fe()[0].dofs_per_vertex +
				 12*this->dof_handler->get_fe()[0].dofs_per_line +
				 6*this->dof_handler->get_fe()[0].dofs_per_quad +
				 this->dof_handler->get_fe()[0].dofs_per_hex),
	  DoFAccessor<hp::DoFHandler<3> >::ExcVectorDoesNotMatch());
  Assert (static_cast<unsigned int>(this->present_level) < this->dof_handler->levels.size(),
          ExcMessage ("DoFHandler not initialized"));

				   // this function really only makes
				   // sense on non-active objects if
				   // all degrees of freedom are
				   // located on vertices, since
				   // otherwise there are degrees of
				   // freedom on sub-objects which are
				   // not allocated for this
				   // non-active thing
  Assert (!this->has_children() ||
	  (this->dof_handler->get_fe()[0].dofs_per_cell ==
	   8*this->dof_handler->get_fe()[0].dofs_per_vertex),
	  DoFAccessor<hp::DoFHandler<3> >::ExcNotActive());
  Assert (static_cast<unsigned int>(this->present_level) < this->dof_handler->levels.size(),
          ExcMessage ("DoFHandler not initialized"));
  
  const unsigned int dofs_per_vertex = this->dof_handler->get_fe()[0].dofs_per_vertex,
		     dofs_per_line   = this->dof_handler->get_fe()[0].dofs_per_line,
		     dofs_per_quad   = this->dof_handler->get_fe()[0].dofs_per_quad,
		     dofs_per_hex    = this->dof_handler->get_fe()[0].dofs_per_hex;
  std::vector<unsigned int>::iterator next = dof_indices.begin();
  for (unsigned int vertex=0; vertex<8; ++vertex)
    for (unsigned int d=0; d<dofs_per_vertex; ++d)
      *next++ = vertex_dof_index(vertex,d);
  for (unsigned int line=0; line<12; ++line)
    for (unsigned int d=0; d<dofs_per_line; ++d)
      *next++ = this->line(line)->dof_index(d);
  for (unsigned int quad=0; quad<6; ++quad)
    for (unsigned int d=0; d<dofs_per_quad; ++d)
      *next++ = this->quad(quad)->dof_index(d);
  for (unsigned int d=0; d<dofs_per_hex; ++d)
    *next++ = dof_index(d);
}




/*------------------------- Functions: DoFCellAccessor -----------------------*/

template <class DH>
inline
DoFCellAccessor<DH>::
DoFCellAccessor (const Triangulation<DH::dimension> *tria,
                 const int                 level,
                 const int                 index,
                 const AccessorData       *local_data)
                :
                DoFObjectAccessor<DH::dimension,DH> (tria,level,index,local_data)
{}


template <class DH>
inline
TriaIterator<DH::dimension,DoFCellAccessor<DH> >
DoFCellAccessor<DH>::neighbor (const unsigned int i) const
{
  TriaIterator<dim,DoFCellAccessor<DH> > q (this->tria,
					    this->neighbor_level (i),
					    this->neighbor_index (i),
					    this->dof_handler);
  
#ifdef DEBUG
  if (q.state() != IteratorState::past_the_end)
    Assert (q->used(), typename TriaAccessor<dim>::ExcUnusedCellAsNeighbor());
#endif
  return q;
}


template <class DH>
inline
TriaIterator<DH::dimension,DoFCellAccessor<DH> >
DoFCellAccessor<DH>::child (const unsigned int i) const
{
  TriaIterator<dim,DoFCellAccessor<DH> > q (this->tria,
					    this->present_level+1,
					    this->child_index (i),
					    this->dof_handler);
  
#ifdef DEBUG
  if (q.state() != IteratorState::past_the_end)
    Assert (q->used(), typename TriaAccessor<dim>::ExcUnusedCellAsChild());
#endif
  return q;
}



template <>
inline
TriaIterator<1, DoFObjectAccessor<0,hp::DoFHandler<1> > >
DoFCellAccessor<hp::DoFHandler<1> >::face (const unsigned int) const
{
  Assert (false, ExcImpossibleInDim(1));
  return TriaIterator<1, DoFObjectAccessor<0,hp::DoFHandler<1> > >();
}



template <>
inline
TriaIterator<2, DoFObjectAccessor<1,hp::DoFHandler<2> > >
DoFCellAccessor<hp::DoFHandler<2> >::face (const unsigned int i) const
{
  return this->line(i);
}



template <>
inline
TriaIterator<3, DoFObjectAccessor<2, hp::DoFHandler<3> > >
DoFCellAccessor<hp::DoFHandler<3> >::face (const unsigned int i) const
{
  return this->quad(i);
}




template <>
inline
void
DoFCellAccessor<DoFHandler<1> >::
get_dof_indices (std::vector<unsigned int> &dof_indices) const
{
  Assert (dof_indices.size() == this->get_fe().dofs_per_cell,
	  ExcVectorDoesNotMatch());

				   // check as in documentation that
				   // cell is either active, or dofs
				   // are only in vertices
  Assert (!this->has_children()
	  ||
	  (this->get_fe().dofs_per_cell ==
	   this->get_fe().dofs_per_vertex * GeometryInfo<1>::vertices_per_cell),
	  ExcMessage ("Cell must either be active, or all DoFs must be in vertices"));
  
  unsigned int *cache = &this->dof_handler->levels[this->present_level]
			->cell_dof_indices_cache[this->present_index * this->get_fe().dofs_per_cell]; 
  for (unsigned int i=0; i<this->get_fe().dofs_per_cell; ++i, ++cache)
    dof_indices[i] = *cache;
}



template <>
inline
void
DoFCellAccessor<hp::DoFHandler<1> >::
get_dof_indices (std::vector<unsigned int> &dof_indices) const
{
				   // no caching for hp::DoFHandler implemented 
  DoFObjectAccessor<dim,hp::DoFHandler<dim> >::get_dof_indices (dof_indices);
}



template <>
inline
void
DoFCellAccessor<DoFHandler<2> >::
get_dof_indices (std::vector<unsigned int> &dof_indices) const
{
  Assert (dof_indices.size() == this->get_fe().dofs_per_cell,
	  ExcVectorDoesNotMatch());

				   // check as in documentation that
				   // cell is either active, or dofs
				   // are only in vertices
  Assert (!this->has_children()
	  ||
	  (this->get_fe().dofs_per_cell ==
	   this->get_fe().dofs_per_vertex * GeometryInfo<2>::vertices_per_cell),
	  ExcMessage ("Cell must either be active, or all DoFs must be in vertices"));
  
  unsigned int *cache = &this->dof_handler->levels[this->present_level]
			->cell_dof_indices_cache[this->present_index * this->get_fe().dofs_per_cell]; 
  for (unsigned int i=0; i<this->get_fe().dofs_per_cell; ++i, ++cache)
    dof_indices[i] = *cache;
}



template <>
inline
void
DoFCellAccessor<hp::DoFHandler<2> >::
get_dof_indices (std::vector<unsigned int> &dof_indices) const
{
				   // no caching for hp::DoFHandler implemented 
  DoFObjectAccessor<dim,hp::DoFHandler<dim> >::get_dof_indices (dof_indices);
}



template <>
inline
void
DoFCellAccessor<DoFHandler<3> >::
get_dof_indices (std::vector<unsigned int> &dof_indices) const
{
  Assert (dof_indices.size() == this->get_fe().dofs_per_cell,
	  ExcVectorDoesNotMatch());

				   // check as in documentation that
				   // cell is either active, or dofs
				   // are only in vertices
  Assert (!this->has_children()
	  ||
	  (this->get_fe().dofs_per_cell ==
	   this->get_fe().dofs_per_vertex * GeometryInfo<3>::vertices_per_cell),
	  ExcMessage ("Cell must either be active, or all DoFs must be in vertices"));
  
  unsigned int *cache = &this->dof_handler->levels[this->present_level]
			->cell_dof_indices_cache[this->present_index * this->get_fe().dofs_per_cell]; 
  for (unsigned int i=0; i<this->get_fe().dofs_per_cell; ++i, ++cache)
    dof_indices[i] = *cache;
}



template <class DH>
inline
const FiniteElement<DH::dimension> &
DoFCellAccessor<DH>::get_fe () const
{
  return this->dof_handler->get_fe();
}



template <class DH>
inline
unsigned int
DoFCellAccessor<DH>::active_fe_index () const
{
  return 0;
}



template <class DH>
inline
void
DoFCellAccessor<DH>::set_active_fe_index (const unsigned int i)
{
  Assert (i == 0, typename DoFAccessor<DH>::ExcInvalidObject());
}





template <>
inline
const FiniteElement<1> &
DoFCellAccessor<hp::DoFHandler<1> >::get_fe () const
{
  return this->dof_handler->get_fe()[active_fe_index ()];
}



template <>
inline
const FiniteElement<2> &
DoFCellAccessor<hp::DoFHandler<2> >::get_fe () const
{
  return this->dof_handler->get_fe()[active_fe_index ()];
}



template <>
inline
const FiniteElement<3> &
DoFCellAccessor<hp::DoFHandler<3> >::get_fe () const
{
  return this->dof_handler->get_fe()[active_fe_index ()];
}



template <>
inline
unsigned int
DoFCellAccessor<hp::DoFHandler<1> >::active_fe_index () const
{
  Assert (static_cast<unsigned int>(this->present_level) < this->dof_handler->levels.size(),
          ExcMessage ("DoFHandler not initialized"));
  Assert (static_cast<std::vector<unsigned int>::size_type>(this->present_index) <
	  this->dof_handler->levels[this->present_level]->active_fe_indices.size (),
	  ExcIndexRange (this->present_index, 0,
			 this->dof_handler->levels[this->present_level]->active_fe_indices.size ()));
  return this->dof_handler->levels[this->present_level]
    ->active_fe_indices[this->present_index];
}



template <>
inline
unsigned int
DoFCellAccessor<hp::DoFHandler<2> >::active_fe_index () const
{
  Assert (static_cast<unsigned int>(this->present_level) < this->dof_handler->levels.size(),
          ExcMessage ("DoFHandler not initialized"));
  Assert (static_cast<std::vector<unsigned int>::size_type>(this->present_index) <
	  this->dof_handler->levels[this->present_level]->active_fe_indices.size (),
	  ExcIndexRange (this->present_index, 0,
			 this->dof_handler->levels[this->present_level]->active_fe_indices.size ()));
  return this->dof_handler->levels[this->present_level]
    ->active_fe_indices[this->present_index];
}



template <>
inline
unsigned int
DoFCellAccessor<hp::DoFHandler<3> >::active_fe_index () const
{
  Assert (static_cast<unsigned int>(this->present_level) < this->dof_handler->levels.size(),
          ExcMessage ("DoFHandler not initialized"));
  Assert (static_cast<std::vector<unsigned int>::size_type>(this->present_index) <
	  this->dof_handler->levels[this->present_level]->active_fe_indices.size (),
	  ExcIndexRange (this->present_index, 0,
			 this->dof_handler->levels[this->present_level]->active_fe_indices.size ()));
  return this->dof_handler->levels[this->present_level]
    ->active_fe_indices[this->present_index];
}



template <>
inline
void
DoFCellAccessor<hp::DoFHandler<1> >::set_active_fe_index (const unsigned int i)
{
  Assert (this->dof_handler != 0,
          DoFAccessor<hp::DoFHandler<1> >::ExcInvalidObject());
  Assert (static_cast<unsigned int>(this->present_level) <
          this->dof_handler->levels.size(),
          ExcMessage ("DoFHandler not initialized"));
  Assert (static_cast<std::vector<unsigned int>::size_type>(this->present_index) <
	  this->dof_handler->levels[this->present_level]->active_fe_indices.size (),
	  ExcIndexRange (this->present_index, 0,
			 this->dof_handler->levels[this->present_level]->active_fe_indices.size ()));
  Assert (i < this->dof_handler->get_fe().size (),
          ExcIndexRange(i, 0, this->dof_handler->get_fe().size ()));
  this->dof_handler->levels[this->present_level]
    ->active_fe_indices[this->present_index] = i;
}



template <>
inline
void
DoFCellAccessor<hp::DoFHandler<2> >::set_active_fe_index (const unsigned int i)
{
  Assert (this->dof_handler != 0,
          DoFAccessor<hp::DoFHandler<2> >::ExcInvalidObject());
  Assert (static_cast<unsigned int>(this->present_level) <
          this->dof_handler->levels.size(),
          ExcMessage ("DoFHandler not initialized"));
  Assert (static_cast<std::vector<unsigned int>::size_type>(this->present_index) <
	  this->dof_handler->levels[this->present_level]->active_fe_indices.size (),
	  ExcIndexRange (this->present_index, 0,
			 this->dof_handler->levels[this->present_level]->active_fe_indices.size ()));
  Assert (i < this->dof_handler->get_fe().size (),
          ExcIndexRange(i, 0, this->dof_handler->get_fe().size ()));
  this->dof_handler->levels[this->present_level]
    ->active_fe_indices[this->present_index] = i;
}



template <>
inline
void
DoFCellAccessor<hp::DoFHandler<3> >::set_active_fe_index (const unsigned int i)
{
  Assert (this->dof_handler != 0,
          DoFAccessor<hp::DoFHandler<3> >::ExcInvalidObject());
  Assert (static_cast<unsigned int>(this->present_level) <
          this->dof_handler->levels.size(),
          ExcMessage ("DoFHandler not initialized"));
  Assert (static_cast<std::vector<unsigned int>::size_type>(this->present_index) <
	  this->dof_handler->levels[this->present_level]->active_fe_indices.size (),
	  ExcIndexRange (this->present_index, 0,
			 this->dof_handler->levels[this->present_level]->active_fe_indices.size ()));
  Assert (i < this->dof_handler->get_fe().size (),
          ExcIndexRange(i, 0, this->dof_handler->get_fe().size ()));
  this->dof_handler->levels[this->present_level]
    ->active_fe_indices[this->present_index] = i;
}



template <>
inline
void
DoFCellAccessor<hp::DoFHandler<3> >::
get_dof_indices (std::vector<unsigned int> &dof_indices) const
{
				   // no caching for hp::DoFHandler implemented 
  DoFObjectAccessor<dim,hp::DoFHandler<dim> >::get_dof_indices (dof_indices);
}



template <class DH>
template <typename number, typename OutputVector>
inline
void
DoFCellAccessor<DH>::
distribute_local_to_global (const Vector<number> &local_source,
			    OutputVector         &global_destination) const
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
  typedef DoFAccessor<DH> BaseClass;
  
  Assert (this->dof_handler != 0,
	  typename BaseClass::ExcInvalidObject());
  Assert (&this->get_fe() != 0,
	  typename BaseClass::ExcInvalidObject());
  Assert (local_source.size() == this->get_fe().dofs_per_cell,
	  typename BaseClass::ExcVectorDoesNotMatch());
  Assert (this->dof_handler->n_dofs() == global_destination.size(),
	  typename BaseClass::ExcVectorDoesNotMatch());

  const unsigned int n_dofs = local_source.size();

//TODO[WB]: This function could me made more efficient. First, it allocates memory, which could be avoided by passing in another argument as a scratch array. second, the elementwise access is really slow if we use PETSc vectors/matrices. This should be fixed eventually
  
				   // get indices of dofs
  std::vector<unsigned int> dofs (n_dofs);
  this->get_dof_indices (dofs);
  
				   // distribute cell vector
  for (unsigned int j=0; j<n_dofs; ++j)
    global_destination(dofs[j]) += local_source(j);
}



template <class DH>
template <typename number, typename OutputMatrix>
inline
void
DoFCellAccessor<DH>::
distribute_local_to_global (const FullMatrix<number> &local_source,
			    OutputMatrix             &global_destination) const
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
  typedef DoFAccessor<DH> BaseClass;
  
  Assert (this->dof_handler != 0,
	  typename BaseClass::ExcInvalidObject());
  Assert (&this->get_fe() != 0,
	  typename BaseClass::ExcInvalidObject());
  Assert (local_source.m() == this->get_fe().dofs_per_cell,
	  typename BaseClass::ExcVectorDoesNotMatch());
  Assert (local_source.m() == local_source.n(),
	  typename BaseClass::ExcMatrixDoesNotMatch());
  Assert (this->dof_handler->n_dofs() == global_destination.m(),
	  typename BaseClass::ExcMatrixDoesNotMatch());
  Assert (global_destination.m() == global_destination.n(),
	  typename BaseClass::ExcMatrixDoesNotMatch());

  const unsigned int n_dofs = local_source.m();

//TODO[WB]: This function could me made more efficient. First, it allocates memory, which could be avoided by passing in another argument as a scratch array. second, the elementwise access is really slow if we use PETSc vectors/matrices. This should be fixed eventually

				   // get indices of dofs
  std::vector<unsigned int> dofs (n_dofs);
  this->get_dof_indices (dofs);
  
				   // distribute cell matrix
  for (unsigned int i=0; i<n_dofs; ++i)
    for (unsigned int j=0; j<n_dofs; ++j)
      global_destination.add(dofs[i], dofs[j], local_source(i,j));
}



#endif
