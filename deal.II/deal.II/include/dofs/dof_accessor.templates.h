//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007 by the deal.II authors
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

DEAL_II_NAMESPACE_OPEN

  

/*------------------------- Functions: DoFAccessor ---------------------------*/


template <int structdim, class DH>
DoFAccessor<structdim,DH>::DoFAccessor ()
                :
		DoFObjectAccessor_Inheritance<structdim, DH::dimension>::BaseClass (0,
										    deal_II_numbers::invalid_unsigned_int,
										    deal_II_numbers::invalid_unsigned_int),
		dof_handler(0)
{
  Assert (false, ExcInvalidObject());
}



template <int structdim, class DH>
inline
DoFAccessor<structdim,DH>::DoFAccessor (const Triangulation<DH::dimension> *tria,
					const int                 level,
					const int                 index,
					const DH                 *dof_handler)
		:
		DoFObjectAccessor_Inheritance<structdim, DH::dimension>::BaseClass (tria,
										    level,
										    index),
                dof_handler(const_cast<DH*>(dof_handler))
{}



template <int structdim, class DH>
inline
void
DoFAccessor<structdim,DH>::set_dof_handler (DH *dh)
{
  Assert (dh != 0, ExcInvalidObject());
  this->dof_handler = dh;
}



template <int structdim, class DH>
inline
const DH &
DoFAccessor<structdim,DH>::get_dof_handler () const
{
  return *this->dof_handler;
}


template <int structdim, class DH>
inline
DoFAccessor<structdim,DH> &
DoFAccessor<structdim,DH>::operator = (const DoFAccessor<structdim,DH> &da)
{
  this->set_dof_handler (da.dof_handler);
  return *this;
}



template <int structdim, class DH>
inline
void
DoFAccessor<structdim,DH>::copy_from (const DoFAccessor<structdim,DH> &a)
{
  BaseClass::copy_from (a);
  set_dof_handler (a.dof_handler);
}



template <int structdim, class DH>
inline
bool
DoFAccessor<structdim,DH>::operator == (const DoFAccessor<structdim,DH> &a) const
{
  Assert (this->dof_handler == a.dof_handler, ExcCantCompareIterators());
  return (BaseClass::operator == (a));
}



template <int structdim, class DH>
inline
bool
DoFAccessor<structdim,DH>::operator != (const DoFAccessor<structdim,DH> &a) const
{
  Assert (this->dof_handler == a.dof_handler, ExcCantCompareIterators());
  return (BaseClass::operator != (a));
}



template <int structdim, class DH>
inline
TriaIterator<DH::dimension,DoFObjectAccessor<structdim,DH> >
DoFAccessor<structdim,DH>::child (const unsigned int i) const
{
  Assert (static_cast<unsigned int>(this->present_level) < this->dof_handler->levels.size(),
          ExcMessage ("DoFHandler not initialized"));

  int next_level;
  if (DH::dimension==structdim)
    next_level = this->present_level+1;
  else
    next_level = 0;
  
  TriaIterator<DH::dimension,DoFObjectAccessor<structdim,DH> > q (this->tria,
								  next_level,
								  this->child_index (i),
								  this->dof_handler);

				   // make sure that we either created
				   // a past-the-end iterator or one
				   // pointing to a used cell
  Assert ((q.state() == IteratorState::past_the_end)
	  ||
	  q->used(),
	  typename TriaAccessor<DH::dimension>::ExcUnusedCellAsChild());

  return q;
}



template <int structdim, class DH>
inline
unsigned int
DoFAccessor<structdim, DH>::vertex_dof_index (const unsigned int vertex,
					      const unsigned int i,
					      const unsigned int fe_index) const
{
  return this->dof_handler->get_vertex_dof_index (this->vertex_index(vertex),
						  fe_index,
						  i);
}



template <int structdim, class DH>
inline
void
DoFAccessor<structdim, DH>::set_vertex_dof_index (const unsigned int vertex,
						  const unsigned int i,
						  const unsigned int index,
						  const unsigned int fe_index) const
{
  this->dof_handler->set_vertex_dof_index (this->vertex_index(vertex),
					   fe_index,
					   i,
					   index);
}



template <int dim, class DH>
inline
unsigned int
DoFAccessor<dim,DH>::dof_index (const unsigned int i,
				const unsigned int fe_index) const
{
				   // access the respective DoF 
  return this->dof_handler
    ->template get_dof_index<dim> (this->present_level,
				   this->present_index,
				   fe_index,
				   i);
}



template <int dim, class DH>
inline
void
DoFAccessor<dim,DH>::set_dof_index (const unsigned int i,
				    const unsigned int index,
				    const unsigned int fe_index) const
{
				   // access the respective DoF
  this->dof_handler
    ->template set_dof_index<dim> (this->present_level,
				   this->present_index,
				   fe_index,
				   i,
				   index);
}



template <int dim, class DH>
inline
unsigned int
DoFAccessor<dim,DH>::n_active_fe_indices () const
{
				   // access the respective DoF
  return this->dof_handler
    ->template n_active_fe_indices<dim> (this->present_level,
					 this->present_index);
}



template <int dim, class DH>
inline
unsigned int
DoFAccessor<dim,DH>::nth_active_fe_index (const unsigned int n) const
{
				   // access the respective DoF
  return this->dof_handler
    ->template nth_active_fe_index<dim> (this->present_level,
					 this->present_index,
					 n);
}



template <int dim, class DH>
inline
bool
DoFAccessor<dim,DH>::fe_index_is_active (const unsigned int fe_index) const
{
				   // access the respective DoF
  return this->dof_handler
    ->template fe_index_is_active<dim> (this->present_level,
					this->present_index,
					fe_index);
}



template <int dim, class DH>
inline
const FiniteElement<DH::dimension> &
DoFAccessor<dim,DH>::get_fe (const unsigned int fe_index) const
{
  Assert (fe_index_is_active (fe_index) == true,
	  ExcMessage ("This function can only be called for active fe indices"));
  
  return this->dof_handler->get_fe()[fe_index];
}



template <>
inline
const FiniteElement<1> &
DoFAccessor<1,dealii::DoFHandler<1> >::get_fe (const unsigned int fe_index) const
{
  Assert (fe_index_is_active (fe_index) == true,
	  ExcMessage ("This function can only be called for active fe indices"));  
  return this->dof_handler->get_fe();
}



template <>
inline
const FiniteElement<2> &
DoFAccessor<1,dealii::DoFHandler<2> >::get_fe (const unsigned int fe_index) const
{
  Assert (fe_index_is_active (fe_index) == true,
	  ExcMessage ("This function can only be called for active fe indices"));  
  return this->dof_handler->get_fe();
}



template <>
inline
const FiniteElement<3> &
DoFAccessor<1,dealii::DoFHandler<3> >::get_fe (const unsigned int fe_index) const
{
  Assert (fe_index_is_active (fe_index) == true,
	  ExcMessage ("This function can only be called for active fe indices"));  
  return this->dof_handler->get_fe();
}



template <>
inline
const FiniteElement<2> &
DoFAccessor<2,dealii::DoFHandler<2> >::get_fe (const unsigned int fe_index) const
{
  Assert (fe_index_is_active (fe_index) == true,
	  ExcMessage ("This function can only be called for active fe indices"));  
  return this->dof_handler->get_fe();
}



template <>
inline
const FiniteElement<3> &
DoFAccessor<2,dealii::DoFHandler<3> >::get_fe (const unsigned int fe_index) const
{
  Assert (fe_index_is_active (fe_index) == true,
	  ExcMessage ("This function can only be called for active fe indices"));  
  return this->dof_handler->get_fe();
}



template <>
inline
const FiniteElement<3> &
DoFAccessor<3,dealii::DoFHandler<3> >::get_fe (const unsigned int fe_index) const
{
  Assert (fe_index_is_active (fe_index) == true,
	  ExcMessage ("This function can only be called for active fe indices"));  
  return this->dof_handler->get_fe();
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
                DoFAccessor<1,DH> (tria, level, index, local_data)
{}




template <class DH>
inline
void
DoFObjectAccessor<1,DH>::get_dof_indices (std::vector<unsigned int> &dof_indices,
					  const unsigned int         fe_index) const
{
  Assert (static_cast<unsigned int>(this->present_level) < this->dof_handler->levels.size(),
          ExcMessage ("DoFHandler not initialized"));

  Assert (this->dof_handler != 0, typename BaseClass::ExcInvalidObject());
  Assert (&this->dof_handler->get_fe() != 0, typename BaseClass::ExcInvalidObject());
  Assert (dof_indices.size() == (2*this->dof_handler->get_fe()[fe_index].dofs_per_vertex +
				 this->dof_handler->get_fe()[fe_index].dofs_per_line),
	  typename BaseClass::ExcVectorDoesNotMatch());

				   // this function really only makes
				   // sense if either a) there are
				   // degrees of freedom defined on
				   // the present object, or b) the
				   // object is non-active objects but
				   // all degrees of freedom are
				   // located on vertices, since
				   // otherwise there are degrees of
				   // freedom on sub-objects which are
				   // not allocated for this
				   // non-active thing
  Assert (this->fe_index_is_active (fe_index)
	  ||
	  (this->dof_handler->get_fe()[fe_index].dofs_per_cell ==
	   2*this->dof_handler->get_fe()[fe_index].dofs_per_vertex),
	  ExcInternalError());
	  
  const unsigned int dofs_per_vertex = this->dof_handler->get_fe()[fe_index].dofs_per_vertex,
		     dofs_per_line   = this->dof_handler->get_fe()[fe_index].dofs_per_line;
  std::vector<unsigned int>::iterator next = dof_indices.begin();
  for (unsigned int vertex=0; vertex<2; ++vertex)
    for (unsigned int d=0; d<dofs_per_vertex; ++d)
      *next++ = this->vertex_dof_index(vertex,d,fe_index);
  for (unsigned int d=0; d<dofs_per_line; ++d)
    *next++ = this->dof_index(d,fe_index);
}





/*------------------------- Functions: DoFObjectAccessor<2,dim> -----------------------*/

template <class DH>
DoFObjectAccessor<2,DH>::
DoFObjectAccessor (const Triangulation<DH::dimension> *tria,
                   const int                 level,
                   const int                 index,
                   const AccessorData       *local_data)
                :
                DoFAccessor<2,DH> (tria, level, index, local_data)
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
      0,
      this->line_index (i),
      this->dof_handler
    );
}



template <class DH>
inline
void
DoFObjectAccessor<2,DH>::get_dof_indices (std::vector<unsigned int> &dof_indices,
					  const unsigned int         fe_index) const
{
  Assert (this->dof_handler != 0,
	  typename BaseClass::ExcInvalidObject());
  Assert (&this->dof_handler->get_fe() != 0,
	  typename BaseClass::ExcInvalidObject());
  Assert (dof_indices.size() == (4*this->dof_handler->get_fe()[fe_index].dofs_per_vertex +
				 4*this->dof_handler->get_fe()[fe_index].dofs_per_line +
				 this->dof_handler->get_fe()[fe_index].dofs_per_quad),
	  typename BaseClass::ExcVectorDoesNotMatch());

				   // this function really only makes
				   // sense if either a) there are
				   // degrees of freedom defined on
				   // the present object, or b) the
				   // object is non-active objects but
				   // all degrees of freedom are
				   // located on vertices, since
				   // otherwise there are degrees of
				   // freedom on sub-objects which are
				   // not allocated for this
				   // non-active thing
  Assert (this->fe_index_is_active (fe_index)
	  ||
	  (this->dof_handler->get_fe()[fe_index].dofs_per_cell ==
	   4*this->dof_handler->get_fe()[fe_index].dofs_per_vertex),
	  ExcInternalError());
	  
  Assert (static_cast<unsigned int>(this->present_level) < this->dof_handler->levels.size(),
          ExcMessage ("DoFHandler not initialized"));
  
  const unsigned int dofs_per_vertex = this->dof_handler->get_fe()[fe_index].dofs_per_vertex,
		     dofs_per_line   = this->dof_handler->get_fe()[fe_index].dofs_per_line,
		     dofs_per_quad   = this->dof_handler->get_fe()[fe_index].dofs_per_quad;
  std::vector<unsigned int>::iterator next = dof_indices.begin();
  for (unsigned int vertex=0; vertex<4; ++vertex)
    for (unsigned int d=0; d<dofs_per_vertex; ++d)
      *next++ = this->vertex_dof_index(vertex,d,fe_index);
				   // now copy dof numbers from the line. for
				   // lines with the wrong orientation (which
				   // might occur in 3d), we have already made
				   // sure that we're ok by picking the correct
				   // vertices (this happens automatically in
				   // the vertex() function). however, if the
				   // line is in wrong orientation, we look at
				   // it in flipped orientation and we will have
				   // to adjust the shape function indices that
				   // we see to correspond to the correct
				   // (face-local) ordering.
  for (unsigned int line=0; line<4; ++line)
    for (unsigned int d=0; d<dofs_per_line; ++d)
      *next++ = this->line(line)->dof_index(this->dof_handler->get_fe()[fe_index].
					    adjust_line_dof_index_for_line_orientation(d,
										       this->line_orientation(line)),fe_index);
  for (unsigned int d=0; d<dofs_per_quad; ++d)
    *next++ = this->dof_index(d,fe_index);
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
                DoFAccessor<3,DH> (tria, level, index, local_data)
{}



template <class DH>
inline
TriaIterator<DH::dimension,DoFObjectAccessor<1,DH> >
DoFObjectAccessor<3,DH>::line (const unsigned int i) const
{
  TriaIterator<dim,TriaObjectAccessor<1,dim> > l = TriaObjectAccessor<3,DH::dimension>::line(i);
  return TriaIterator<dim,DoFObjectAccessor<1,DH> >
    (
      this->tria,
      0,
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
      0,
      this->quad_index (i),
      this->dof_handler
    );
}



template <class DH>
inline
void
DoFObjectAccessor<3,DH>::get_dof_indices (std::vector<unsigned int> &dof_indices,
					  const unsigned int         fe_index) const
{
  Assert (this->dof_handler != 0,
	  typename BaseClass::ExcInvalidObject());
  Assert (&this->dof_handler->get_fe() != 0,
	  typename BaseClass::ExcInvalidObject());
  Assert (dof_indices.size() == (8*this->dof_handler->get_fe()[fe_index].dofs_per_vertex +
				 12*this->dof_handler->get_fe()[fe_index].dofs_per_line +
				 6*this->dof_handler->get_fe()[fe_index].dofs_per_quad +
				 this->dof_handler->get_fe()[fe_index].dofs_per_hex),
	  typename BaseClass::ExcVectorDoesNotMatch());
  Assert (static_cast<unsigned int>(this->present_level) < this->dof_handler->levels.size(),
          ExcMessage ("DoFHandler not initialized"));

				   // this function really only makes
				   // sense if either a) there are
				   // degrees of freedom defined on
				   // the present object, or b) the
				   // object is non-active objects but
				   // all degrees of freedom are
				   // located on vertices, since
				   // otherwise there are degrees of
				   // freedom on sub-objects which are
				   // not allocated for this
				   // non-active thing
  Assert (this->fe_index_is_active (fe_index)
	  ||
	  (this->dof_handler->get_fe()[fe_index].dofs_per_cell ==
	   8*this->dof_handler->get_fe()[fe_index].dofs_per_vertex),
	  ExcInternalError());

  Assert (static_cast<unsigned int>(this->present_level) < this->dof_handler->levels.size(),
          ExcMessage ("DoFHandler not initialized"));
  
  const unsigned int dofs_per_vertex = this->dof_handler->get_fe()[fe_index].dofs_per_vertex,
		     dofs_per_line   = this->dof_handler->get_fe()[fe_index].dofs_per_line,
		     dofs_per_quad   = this->dof_handler->get_fe()[fe_index].dofs_per_quad,
		     dofs_per_hex    = this->dof_handler->get_fe()[fe_index].dofs_per_hex;
  std::vector<unsigned int>::iterator next = dof_indices.begin();
  for (unsigned int vertex=0; vertex<8; ++vertex)
    for (unsigned int d=0; d<dofs_per_vertex; ++d)
      *next++ = this->vertex_dof_index(vertex,d,fe_index);
				   // now copy dof numbers from the line. for
				   // lines with the wrong orientation, we have
				   // already made sure that we're ok by picking
				   // the correct vertices (this happens
				   // automatically in the vertex()
				   // function). however, if the line is in
				   // wrong orientation, we look at it in
				   // flipped orientation and we will have to
				   // adjust the shape function indices that we
				   // see to correspond to the correct
				   // (cell-local) ordering.
  for (unsigned int line=0; line<12; ++line)
    for (unsigned int d=0; d<dofs_per_line; ++d)
      *next++ = this->line(line)->dof_index(this->dof_handler->get_fe()[fe_index].
					    adjust_line_dof_index_for_line_orientation(d,
										       this->line_orientation(line)),fe_index);
				   // now copy dof numbers from the face. for
				   // faces with the wrong orientation, we
				   // have already made sure that we're ok by
				   // picking the correct lines and vertices
				   // (this happens automatically in the
				   // line() and vertex() functions). however,
				   // if the face is in wrong orientation, we
				   // look at it in flipped orientation and we
				   // will have to adjust the shape function
				   // indices that we see to correspond to the
				   // correct (cell-local) ordering. The same
				   // applies, if the face_rotation or
				   // face_orientation is non-standard
  for (unsigned int quad=0; quad<6; ++quad)
    for (unsigned int d=0; d<dofs_per_quad; ++d)
      *next++ = this->quad(quad)->dof_index(this->dof_handler->get_fe()[fe_index].
					      adjust_quad_dof_index_for_face_orientation(d,
											 this->face_orientation(quad),
											 this->face_flip(quad),
											 this->face_rotation(quad)),fe_index);
  for (unsigned int d=0; d<dofs_per_hex; ++d)
    *next++ = this->dof_index(d,fe_index);
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
  DoFObjectAccessor<dim,hp::DoFHandler<dim> >::get_dof_indices (dof_indices,
								this->active_fe_index());
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
  DoFObjectAccessor<dim,hp::DoFHandler<dim> >::get_dof_indices (dof_indices,
								this->active_fe_index());
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



template <>
inline
void
DoFCellAccessor<hp::DoFHandler<3> >::
get_dof_indices (std::vector<unsigned int> &dof_indices) const
{
				   // no caching for hp::DoFHandler implemented 
  DoFObjectAccessor<dim,hp::DoFHandler<dim> >::get_dof_indices (dof_indices,
								this->active_fe_index());
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
  Assert (i == 0, typename BaseClass::BaseClass::ExcInvalidObject());
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
          BaseClass::ExcInvalidObject());
  Assert (static_cast<unsigned int>(this->present_level) <
          this->dof_handler->levels.size(),
          ExcMessage ("DoFHandler not initialized"));
  Assert (static_cast<std::vector<unsigned int>::size_type>(this->present_index) <
	  this->dof_handler->levels[this->present_level]->active_fe_indices.size (),
	  ExcIndexRange (this->present_index, 0,
			 this->dof_handler->levels[this->present_level]->active_fe_indices.size ()));
  this->dof_handler->levels[this->present_level]
    ->active_fe_indices[this->present_index] = i;
}



template <>
inline
void
DoFCellAccessor<hp::DoFHandler<2> >::set_active_fe_index (const unsigned int i)
{
  Assert (this->dof_handler != 0,
          BaseClass::ExcInvalidObject());
  Assert (static_cast<unsigned int>(this->present_level) <
          this->dof_handler->levels.size(),
          ExcMessage ("DoFHandler not initialized"));
  Assert (static_cast<std::vector<unsigned int>::size_type>(this->present_index) <
	  this->dof_handler->levels[this->present_level]->active_fe_indices.size (),
	  ExcIndexRange (this->present_index, 0,
			 this->dof_handler->levels[this->present_level]->active_fe_indices.size ()));
  this->dof_handler->levels[this->present_level]
    ->active_fe_indices[this->present_index] = i;
}



template <>
inline
void
DoFCellAccessor<hp::DoFHandler<3> >::set_active_fe_index (const unsigned int i)
{
  Assert (this->dof_handler != 0,
          BaseClass::ExcInvalidObject());
  Assert (static_cast<unsigned int>(this->present_level) <
          this->dof_handler->levels.size(),
          ExcMessage ("DoFHandler not initialized"));
  Assert (static_cast<std::vector<unsigned int>::size_type>(this->present_index) <
	  this->dof_handler->levels[this->present_level]->active_fe_indices.size (),
	  ExcIndexRange (this->present_index, 0,
			 this->dof_handler->levels[this->present_level]->active_fe_indices.size ()));
  this->dof_handler->levels[this->present_level]
    ->active_fe_indices[this->present_index] = i;
}



template <class DH>
template <typename number, typename OutputVector>
inline
void
DoFCellAccessor<DH>::
distribute_local_to_global (const Vector<number> &local_source,
			    OutputVector         &global_destination) const
{
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


DEAL_II_NAMESPACE_CLOSE

#endif
