//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008 by the deal.II authors
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
#include <hp/dof_levels.h>
#include <grid/tria_iterator.h>
#include <grid/tria_iterator.templates.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN


/*------------------------- Functions: DoFAccessor ---------------------------*/


template <int structdim, class DH>
DoFAccessor<structdim,DH>::
DoFAccessor ()
{
  Assert (false, ExcInvalidObject());
}



template <int structdim, class DH>
inline
DoFAccessor<structdim,DH>::
DoFAccessor (const Triangulation<DH::dimension,DH::space_dimension> *tria,
	     const int                 level,
	     const int                 index,
	     const DH                 *dof_handler)
		:
		internal::DoFAccessor::Inheritance<structdim,DH::dimension,DH::space_dimension>::BaseClass (tria,
													    level,
													    index),
                dof_handler(const_cast<DH*>(dof_handler))
{}



template <int structdim, class DH>
template <int structdim2, int dim2, int spacedim2>
DoFAccessor<structdim,DH>::
DoFAccessor (const InvalidAccessor<structdim2,dim2,spacedim2> &)
{
  Assert (false, ExcInvalidObject());
}



template <int structdim, class DH>
template <int dim2, class DH2>
DoFAccessor<structdim,DH>::
DoFAccessor (const DoFAccessor<dim2, DH2> &)
{
  Assert (false, ExcInvalidObject());
}



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
TriaIterator<DoFAccessor<structdim,DH> >
DoFAccessor<structdim,DH>::child (const unsigned int i) const
{
  Assert (static_cast<unsigned int>(this->level()) < this->dof_handler->levels.size(),
          ExcMessage ("DoFHandler not initialized"));

  int next_level;
  if (DH::dimension==structdim)
    next_level = this->level()+1;
  else
    next_level = 0;
  
  TriaIterator<DoFAccessor<structdim,DH> > q (this->tria,
					      next_level,
					      this->child_index (i),
					      this->dof_handler);
  
				   // make sure that we either created
				   // a past-the-end iterator or one
				   // pointing to a used cell
  Assert ((q.state() == IteratorState::past_the_end)
	  ||
	  q->used(),
	  TriaAccessorExceptions::ExcUnusedCellAsChild());

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
    ->template get_dof_index<dim> (this->level(),
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
    ->template set_dof_index<dim> (this->level(),
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
    ->template n_active_fe_indices<dim> (this->level(),
					 this->present_index);
}



template <int dim, class DH>
inline
unsigned int
DoFAccessor<dim,DH>::nth_active_fe_index (const unsigned int n) const
{
				   // access the respective DoF
  return this->dof_handler
    ->template nth_active_fe_index<dim> (this->level(),
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
    ->template fe_index_is_active<dim> (this->level(),
					this->present_index,
					fe_index);
}



namespace internal
{
  namespace DoFAccessor
  {
    template <int dim, int spacedim>
    inline
    const FiniteElement<dim,spacedim> &
    get_fe (const FiniteElement<dim,spacedim> &fe,
	    const unsigned int)
    {
      return fe;
    }


  
    template <int dim, int spacedim>
    inline
    const FiniteElement<dim,spacedim> &
    get_fe (const dealii::hp::FECollection<dim,spacedim> &fe,
	    const unsigned int                            index)
    {
      return fe[index];
    }
  }
}


template <int dim, class DH>
inline
const FiniteElement<DH::dimension,DH::space_dimension> &
DoFAccessor<dim,DH>::get_fe (const unsigned int fe_index) const
{
  Assert (fe_index_is_active (fe_index) == true,
	  ExcMessage ("This function can only be called for active fe indices"));
  
  return internal::DoFAccessor::get_fe (this->dof_handler->get_fe(), fe_index);
}



namespace internal
{
  namespace DoFAccessor
  {
    template <class DH>
    void get_dof_indices (const dealii::DoFAccessor<1,DH>   &accessor,
			  std::vector<unsigned int> &dof_indices,
			  const unsigned int         fe_index)
    {
      const unsigned int dofs_per_vertex = accessor.get_fe(fe_index).dofs_per_vertex,
			 dofs_per_line   = accessor.get_fe(fe_index).dofs_per_line;
      std::vector<unsigned int>::iterator next = dof_indices.begin();
      for (unsigned int vertex=0; vertex<2; ++vertex)
	for (unsigned int d=0; d<dofs_per_vertex; ++d)
	  *next++ = accessor.vertex_dof_index(vertex,d,fe_index);
      for (unsigned int d=0; d<dofs_per_line; ++d)
	*next++ = accessor.dof_index(d,fe_index);
    }


    
    template <class DH>
    void get_dof_indices (const dealii::DoFAccessor<2,DH>   &accessor,
			  std::vector<unsigned int> &dof_indices,
			  const unsigned int         fe_index)
    {
      const unsigned int dofs_per_vertex = accessor.get_fe(fe_index).dofs_per_vertex,
			 dofs_per_line   = accessor.get_fe(fe_index).dofs_per_line,
			 dofs_per_quad   = accessor.get_fe(fe_index).dofs_per_quad;
      std::vector<unsigned int>::iterator next = dof_indices.begin();
      for (unsigned int vertex=0; vertex<4; ++vertex)
	for (unsigned int d=0; d<dofs_per_vertex; ++d)
	  *next++ = accessor.vertex_dof_index(vertex,d,fe_index);
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
	  *next++ = accessor.line(line)->dof_index(accessor.get_fe(fe_index).
						   adjust_line_dof_index_for_line_orientation(d,
											      accessor.line_orientation(line)),
						   fe_index);
      for (unsigned int d=0; d<dofs_per_quad; ++d)
	*next++ = accessor.dof_index(d,fe_index);
    }

    

    template <class DH>
    void get_dof_indices (const dealii::DoFAccessor<3,DH>   &accessor,
			  std::vector<unsigned int> &dof_indices,
			  const unsigned int         fe_index)
    {
      const unsigned int dofs_per_vertex = accessor.get_fe(fe_index).dofs_per_vertex,
			 dofs_per_line   = accessor.get_fe(fe_index).dofs_per_line,
			 dofs_per_quad   = accessor.get_fe(fe_index).dofs_per_quad,
			 dofs_per_hex    = accessor.get_fe(fe_index).dofs_per_hex;
      std::vector<unsigned int>::iterator next = dof_indices.begin();
      for (unsigned int vertex=0; vertex<8; ++vertex)
	for (unsigned int d=0; d<dofs_per_vertex; ++d)
	  *next++ = accessor.vertex_dof_index(vertex,d,fe_index);
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
	  *next++ = accessor.line(line)->dof_index(accessor.get_fe(fe_index).
						   adjust_line_dof_index_for_line_orientation(d,
											      accessor.line_orientation(line)),fe_index);
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
	  *next++ = accessor.quad(quad)->dof_index(accessor.get_fe(fe_index).
						   adjust_quad_dof_index_for_face_orientation(d,
											      accessor.face_orientation(quad),
											      accessor.face_flip(quad),
											      accessor.face_rotation(quad)),
						   fe_index);
      for (unsigned int d=0; d<dofs_per_hex; ++d)
	*next++ = accessor.dof_index(d,fe_index);
    }
  }
}


template <int structdim, class DH>
inline
void
DoFAccessor<structdim,DH>::get_dof_indices (std::vector<unsigned int> &dof_indices,
					    const unsigned int         fe_index) const
{
  Assert (static_cast<unsigned int>(this->level()) < this->dof_handler->levels.size(),
          ExcMessage ("DoFHandler not initialized"));

  Assert (this->dof_handler != 0, ExcInvalidObject());
  Assert (&this->dof_handler->get_fe() != 0, ExcInvalidObject());
  Assert (static_cast<unsigned int>(this->level()) < this->dof_handler->levels.size(),
          ExcMessage ("DoFHandler not initialized"));

  switch (structdim)
    {
      case 1:
	    Assert (dof_indices.size() ==
		    (2*this->dof_handler->get_fe()[fe_index].dofs_per_vertex +
		     this->dof_handler->get_fe()[fe_index].dofs_per_line),
		    ExcVectorDoesNotMatch());
	    break;
      case 2:
	    Assert (dof_indices.size() ==
		    (4*this->dof_handler->get_fe()[fe_index].dofs_per_vertex +
		     4*this->dof_handler->get_fe()[fe_index].dofs_per_line +
		     this->dof_handler->get_fe()[fe_index].dofs_per_quad),
		    ExcVectorDoesNotMatch());
	    break;
      case 3:
	    Assert (dof_indices.size() ==
		    (8*this->dof_handler->get_fe()[fe_index].dofs_per_vertex +
		     12*this->dof_handler->get_fe()[fe_index].dofs_per_line +
		     6*this->dof_handler->get_fe()[fe_index].dofs_per_quad +
		     this->dof_handler->get_fe()[fe_index].dofs_per_hex),
		    ExcVectorDoesNotMatch());
	    break;
      default:
	    Assert (false, ExcNotImplemented());
    }
  

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
	   GeometryInfo<structdim>::vertices_per_cell *
	   this->dof_handler->get_fe()[fe_index].dofs_per_vertex),
	  ExcInternalError());

				   // now do the actual work
  internal::DoFAccessor::get_dof_indices (*this, dof_indices, fe_index);
}





template <int structdim, class DH>
inline
typename internal::DoFHandler::Iterators<DH>::line_iterator
DoFAccessor<structdim,DH>::line (const unsigned int i) const
{
  Assert (structdim > 1, ExcImpossibleInDim(structdim));
				   // checking of 'i' happens in
				   // line_index(i)
  
  return typename internal::DoFHandler::Iterators<DH>::line_iterator
    (
      this->tria,
      0,  // only sub-objects are allowed, which have no level
      this->line_index(i),
      this->dof_handler
    );
}


template <int structdim, class DH>
inline
typename internal::DoFHandler::Iterators<DH>::quad_iterator
DoFAccessor<structdim,DH>::quad (const unsigned int i) const
{
  Assert (structdim > 2, ExcImpossibleInDim(structdim));
				   // checking of 'i' happens in
				   // quad_index(i)

  return typename internal::DoFHandler::Iterators<DH>::quad_iterator
    (
      this->tria,
      0,  // only sub-objects are allowed, which have no level
      this->quad_index (i),
      this->dof_handler
    );
}





/*------------------------- Functions: DoFCellAccessor -----------------------*/


namespace internal
{
  namespace DoFCellAccessor
  {
				     // make sure we refer to class
				     // dealii::DoFCellAccessor, not
				     // namespace
				     // internal::DoFCellAccessor
    using dealii::DoFCellAccessor;
    using dealii::DoFHandler;
    
/**
 * A class with the same purpose as the similarly named class of the
 * Triangulation class. See there for more information.
 */
    struct Implementation
    {
					 /**
					  * Implement the updating of the
					  * cache. Currently not
					  * implemented for hp::DoFHandler
					  * objects.
					  */
	template <int spacedim>
	static
	void
	update_cell_dof_indices_cache (const DoFCellAccessor<DoFHandler<1,spacedim> > &accessor)
	  {
					     // check as in documentation that
					     // cell is either active, or dofs
					     // are only in vertices. otherwise
					     // simply don't update the cache at
					     // all. the get_dof_indices
					     // function will then make sure we
					     // don't access the invalid data
	    if (accessor.has_children()
		&&
		(accessor.get_fe().dofs_per_cell !=
		 accessor.get_fe().dofs_per_vertex * GeometryInfo<1>::vertices_per_cell))
	      return;
	
	    const unsigned int dofs_per_vertex = accessor.get_fe().dofs_per_vertex,
			       dofs_per_line   = accessor.get_fe().dofs_per_line,
			       dofs_per_cell   = accessor.get_fe().dofs_per_cell;
	
					     // make sure the cache is at least
					     // as big as we need it when
					     // writing to the last element of
					     // this cell
	    Assert (accessor.present_index * dofs_per_cell + dofs_per_cell
		    <=
		    accessor.dof_handler->levels[accessor.present_level]
		    ->cell_dof_indices_cache.size(),
		    ExcInternalError());
	
	    std::vector<unsigned int>::iterator next
	      = (accessor.dof_handler->levels[accessor.present_level]
		 ->cell_dof_indices_cache.begin() + accessor.present_index * dofs_per_cell);
	
	    for (unsigned int vertex=0; vertex<2; ++vertex)
	      for (unsigned int d=0; d<dofs_per_vertex; ++d)
		*next++ = accessor.vertex_dof_index(vertex,d);
	    for (unsigned int d=0; d<dofs_per_line; ++d)
	      *next++ = accessor.dof_index(d);
	  }



	template <int spacedim>
	static
	void
	update_cell_dof_indices_cache (const DoFCellAccessor<DoFHandler<2,spacedim> > &accessor)
	  {
					     // check as in documentation that
					     // cell is either active, or dofs
					     // are only in vertices. otherwise
					     // simply don't update the cache at
					     // all. the get_dof_indices
					     // function will then make sure we
					     // don't access the invalid data
	    if (accessor.has_children()
		&&
		(accessor.get_fe().dofs_per_cell !=
		 accessor.get_fe().dofs_per_vertex * GeometryInfo<2>::vertices_per_cell))
	      return;
  
	    const unsigned int dofs_per_vertex = accessor.get_fe().dofs_per_vertex,
			       dofs_per_line   = accessor.get_fe().dofs_per_line,
			       dofs_per_quad   = accessor.get_fe().dofs_per_quad,
			       dofs_per_cell   = accessor.get_fe().dofs_per_cell;

					     // make sure the cache is at least
					     // as big as we need it when
					     // writing to the last element of
					     // this cell
	    Assert (accessor.present_index * dofs_per_cell + dofs_per_cell
		    <=
		    accessor.dof_handler->levels[accessor.present_level]
		    ->cell_dof_indices_cache.size(),
		    ExcInternalError());

	    std::vector<unsigned int>::iterator next
	      = (accessor.dof_handler->levels[accessor.present_level]
		 ->cell_dof_indices_cache.begin() + accessor.present_index * dofs_per_cell);

	    for (unsigned int vertex=0; vertex<4; ++vertex)
	      for (unsigned int d=0; d<dofs_per_vertex; ++d)
		*next++ = accessor.vertex_dof_index(vertex,d);
	    for (unsigned int line=0; line<4; ++line)
	      for (unsigned int d=0; d<dofs_per_line; ++d)
		*next++ = accessor.line(line)->dof_index(d);
	    for (unsigned int d=0; d<dofs_per_quad; ++d)
	      *next++ = accessor.dof_index(d);
	  }



	template <int spacedim>
	static
	void
	update_cell_dof_indices_cache (const DoFCellAccessor<DoFHandler<3,spacedim> > &accessor)
	  {
					     // check as in documentation that
					     // cell is either active, or dofs
					     // are only in vertices. otherwise
					     // simply don't update the cache at
					     // all. the get_dof_indices
					     // function will then make sure we
					     // don't access the invalid data
	    if (accessor.has_children()
		&&
		(accessor.get_fe().dofs_per_cell !=
		 accessor.get_fe().dofs_per_vertex * GeometryInfo<3>::vertices_per_cell))
	      return;
  
	    const unsigned int dofs_per_vertex = accessor.get_fe().dofs_per_vertex,
			       dofs_per_line   = accessor.get_fe().dofs_per_line,
			       dofs_per_quad   = accessor.get_fe().dofs_per_quad,
			       dofs_per_hex    = accessor.get_fe().dofs_per_hex,
			       dofs_per_cell   = accessor.get_fe().dofs_per_cell;

					     // make sure the cache is at least
					     // as big as we need it when
					     // writing to the last element of
					     // this cell
	    Assert (accessor.present_index * dofs_per_cell + dofs_per_cell
		    <=
		    accessor.dof_handler->levels[accessor.present_level]
		    ->cell_dof_indices_cache.size(),
		    ExcInternalError());

	    std::vector<unsigned int>::iterator next
	      = (accessor.dof_handler->levels[accessor.present_level]
		 ->cell_dof_indices_cache.begin() + accessor.present_index * dofs_per_cell);

	    for (unsigned int vertex=0; vertex<8; ++vertex)
	      for (unsigned int d=0; d<dofs_per_vertex; ++d)
		*next++ = accessor.vertex_dof_index(vertex,d);
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
		*next++ = accessor.line(line)->dof_index(accessor.dof_handler->get_fe().
							 adjust_line_dof_index_for_line_orientation(d,
												    accessor.line_orientation(line)));
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
		*next++ = accessor.quad(quad)->dof_index(accessor.dof_handler->get_fe().
							 adjust_quad_dof_index_for_face_orientation(d,
												    accessor.face_orientation(quad),
												    accessor.face_flip(quad),
												    accessor.face_rotation(quad)));
	    for (unsigned int d=0; d<dofs_per_hex; ++d)
	      *next++ = accessor.dof_index(d);
	  }


					 // implementation for the case of
					 // hp::DoFHandler objects. it's
					 // not implemented there, for no
					 // space dimension
	template <int dim, int spacedim>
	static
	void
	update_cell_dof_indices_cache (const DoFCellAccessor<dealii::hp::DoFHandler<dim,spacedim> > &)
	  {
//TODO[WB]: should implement a dof indices cache for hp as well
  
					     // not implemented, but should also
					     // not be called
	    Assert (false, ExcNotImplemented());
	  }


					 /**
					  * A function that collects the
					  * global indices of degrees of
					  * freedom. This function works
					  * for ::DoFHandler and all
					  * template arguments and copies
					  * the data out of the cache that
					  * we hold for each cell.
					  */
	template <int dim, int spacedim>
	static
	void
	get_dof_indices (const DoFCellAccessor<DoFHandler<dim,spacedim> > &accessor,
			 std::vector<unsigned int>                        &dof_indices)
	  {
	    typedef
	      dealii::DoFAccessor<dim,DoFHandler<dim,spacedim> >
	      BaseClass;
	    Assert (dof_indices.size() == accessor.get_fe().dofs_per_cell,
		    typename BaseClass::ExcVectorDoesNotMatch());

					     // check as in documentation that
					     // cell is either active, or dofs
					     // are only in vertices
	    Assert (!accessor.has_children()
		    ||
		    (accessor.get_fe().dofs_per_cell ==
		     accessor.get_fe().dofs_per_vertex * GeometryInfo<dim>::vertices_per_cell),
		    ExcMessage ("Cell must either be active, or all DoFs must be in vertices"));
  
	    unsigned int *cache = &accessor.dof_handler->levels[accessor.level()]
				  ->cell_dof_indices_cache[accessor.present_index * accessor.get_fe().dofs_per_cell]; 
	    for (unsigned int i=0; i<accessor.get_fe().dofs_per_cell; ++i, ++cache)
	      dof_indices[i] = *cache;
	  }      

					 /**
					  * Same function as above except
					  * that it works for
					  * hp::DoFHandler objects that do
					  * not have a cache for the local
					  * DoF indices.
					  */
	template <int dim, int spacedim>
	static
	void
	get_dof_indices (const DoFCellAccessor<dealii::hp::DoFHandler<dim,spacedim> > &accessor,
			 std::vector<unsigned int>                            &dof_indices)
	  {
					     // no caching for hp::DoFHandler implemented 
	    accessor.dealii::DoFAccessor<dim,dealii::hp::DoFHandler<dim,spacedim> >::get_dof_indices (dof_indices,
										      accessor.active_fe_index());
	  }


					 /**
					  * Do what the active_fe_index
					  * function in the parent class
					  * is supposed to do.
					  */
	template <int dim, int spacedim>
	static
	unsigned int
	active_fe_index (const DoFCellAccessor<DoFHandler<dim,spacedim> > &)
	  {
					     // ::DoFHandler only supports a
					     // single active fe with index
					     // zero
	    return 0;
	  }



	template <int dim, int spacedim>
	static
	unsigned int
	active_fe_index (const DoFCellAccessor<dealii::hp::DoFHandler<dim,spacedim> > &accessor)
	  {
	    Assert (static_cast<unsigned int>(accessor.level()) < accessor.dof_handler->levels.size(),
		    ExcMessage ("DoFHandler not initialized"));
	    Assert (static_cast<std::vector<unsigned int>::size_type>(accessor.present_index) <
		    accessor.dof_handler->levels[accessor.level()]->active_fe_indices.size (),
		    ExcIndexRange (accessor.present_index, 0,
				   accessor.dof_handler->levels[accessor.level()]->active_fe_indices.size ()));
	    return accessor.dof_handler->levels[accessor.level()]
	      ->active_fe_indices[accessor.present_index];
	  }



					 /**
					  * Do what the
					  * set_active_fe_index function
					  * in the parent class is
					  * supposed to do.
					  */
	template <int dim, int spacedim>
	static
	void
	set_active_fe_index (const DoFCellAccessor<DoFHandler<dim,spacedim> > &,
			     const unsigned int                                i)
	  {
					     // ::DoFHandler only supports a
					     // single active fe with index
					     // zero
	    typedef
	      dealii::DoFAccessor<dim,DoFHandler<dim,spacedim> >
	      BaseClass;
	    Assert (i == 0, typename BaseClass::ExcInvalidObject());
	  }



	template <int dim, int spacedim>
	static
	void
	set_active_fe_index (DoFCellAccessor<dealii::hp::DoFHandler<dim,spacedim> > &accessor,
			     const unsigned int                              i)
	  {
	    typedef
	      dealii::DoFAccessor<dim,DoFHandler<dim,spacedim> >
	      BaseClass;
	    Assert (accessor.dof_handler != 0,
		    typename BaseClass::ExcInvalidObject());
	    Assert (static_cast<unsigned int>(accessor.level()) <
		    accessor.dof_handler->levels.size(),
		    ExcMessage ("DoFHandler not initialized"));
	    Assert (static_cast<std::vector<unsigned int>::size_type>(accessor.present_index) <
		    accessor.dof_handler->levels[accessor.level()]->active_fe_indices.size (),
		    ExcIndexRange (accessor.present_index, 0,
				   accessor.dof_handler->levels[accessor.level()]->active_fe_indices.size ()));
	    accessor.dof_handler->levels[accessor.level()]
	      ->active_fe_indices[accessor.present_index] = i;
	  }
    };
  }
}




template <class DH>
inline
DoFCellAccessor<DH>::
DoFCellAccessor (const Triangulation<DH::dimension,DH::space_dimension> *tria,
                 const int                 level,
                 const int                 index,
                 const AccessorData       *local_data)
                :
                DoFAccessor<DH::dimension,DH> (tria,level,index,local_data)
{}



template <class DH>
template <int structdim2, int dim2, int spacedim2>
inline
DoFCellAccessor<DH>::
DoFCellAccessor (const InvalidAccessor<structdim2,dim2,spacedim2> &)
{
  Assert (false, typename BaseClass::ExcInvalidObject());
}



template <class DH>
template <int dim2, class DH2>
inline
DoFCellAccessor<DH>::
DoFCellAccessor (const DoFAccessor<dim2,DH2> &)
{
  Assert (false, typename BaseClass::ExcInvalidObject());
}


template <class DH>
inline
typename internal::DoFHandler::Iterators<DH>::cell_iterator
DoFCellAccessor<DH>::neighbor (const unsigned int i) const
{
  typename internal::DoFHandler::Iterators<DH>::cell_iterator
    q (this->tria,
       this->neighbor_level (i),
       this->neighbor_index (i),
       this->dof_handler);
  
#ifdef DEBUG
  if (q.state() != IteratorState::past_the_end)
    Assert (q->used(), TriaAccessorExceptions::ExcUnusedCellAsNeighbor());
#endif
  return q;
}


template <class DH>
inline
typename internal::DoFHandler::Iterators<DH>::cell_iterator
DoFCellAccessor<DH>::child (const unsigned int i) const
{
  typename internal::DoFHandler::Iterators<DH>::cell_iterator
    q (this->tria,
       this->level()+1,
       this->child_index (i),
       this->dof_handler);
  
#ifdef DEBUG
  if (q.state() != IteratorState::past_the_end)
    Assert (q->used(), TriaAccessorExceptions::ExcUnusedCellAsChild());
#endif
  return q;
}



template <class DH>
typename internal::DoFHandler::Iterators<DH>::face_iterator
DoFCellAccessor<DH>::face (const unsigned int i) const
{
  Assert (i<GeometryInfo<dim>::faces_per_cell, ExcIndexRange (i, 0, GeometryInfo<dim>::faces_per_cell));
  Assert (static_cast<unsigned int>(this->level()) < this->dof_handler->levels.size(),
          ExcMessage ("DoFHandler not initialized"));

  const unsigned int dim = DH::dimension;  
  Assert (dim > 1, ExcImpossibleInDim(1));

  switch (dim)
    {
      case 2:
	    return typename internal::DoFHandler::Iterators<DH>::face_iterator
	      (this->tria,
	       0,
	       this->line_index (i),
	       this->dof_handler);
	    
      case 3:
	    return typename internal::DoFHandler::Iterators<DH>::face_iterator
	      (this->tria,
	       0,
	       this->quad_index (i),
	       this->dof_handler);
	    
      default:
	    Assert (false, ExcNotImplemented());
	    return typename internal::DoFHandler::Iterators<DH>::face_iterator();
    }
}



template <class DH>
inline
void
DoFCellAccessor<DH>::
get_dof_indices (std::vector<unsigned int> &dof_indices) const
{
  internal::DoFCellAccessor::Implementation::get_dof_indices (*this, dof_indices);
}



template <class DH>
inline
const FiniteElement<DH::dimension,DH::space_dimension> &
DoFCellAccessor<DH>::get_fe () const
{
  return internal::DoFAccessor::get_fe (this->dof_handler->get_fe(), active_fe_index());
}



template <class DH>
inline
unsigned int
DoFCellAccessor<DH>::active_fe_index () const
{
  return internal::DoFCellAccessor::Implementation::active_fe_index (*this);
}



template <class DH>
inline
void
DoFCellAccessor<DH>::set_active_fe_index (const unsigned int i)
{
  internal::DoFCellAccessor::Implementation::set_active_fe_index (*this, i);
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
