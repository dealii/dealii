//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <lac/vector.h>
#include <lac/block_vector.h>
#include <lac/petsc_vector.h>
#include <lac/petsc_block_vector.h>
#include <lac/sparse_matrix.h>

#include <dofs/dof_accessor.h>
#include <dofs/dof_accessor.templates.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_levels.h>
#include <dofs/hp_dof_handler.h>
#include <grid/tria_boundary.h>
#include <grid/tria_iterator.h>
#include <grid/tria_iterator.templates.h>
#include <fe/fe.h>

#include <vector>


/*------------------------- Functions: DoFObjectAccessor<1,dim> -----------------------*/


template <class DH>
void DoFObjectAccessor<1, DH>::set_dof_index (const unsigned int i,
                                              const unsigned int index,
					      const unsigned int fe_index) const
{
  Assert (fe_index == DoFHandler<1>::default_fe_index,
	  ExcMessage ("Only the default FE index is allowed for non-hp DoFHandler objects"));
  typedef DoFAccessor<DH> BaseClass;

  Assert (this->dof_handler != 0,
	  typename BaseClass::ExcInvalidObject());
				   // make sure a FE has been selected
				   // and enough room was reserved
  Assert (&this->dof_handler->get_fe() != 0,
	  typename BaseClass::ExcInvalidObject());
  Assert (i<this->dof_handler->get_fe().dofs_per_line,
	  ExcIndexRange (i, 0, this->dof_handler->get_fe().dofs_per_line));

  this->dof_handler->levels[this->present_level]
    ->line_dofs[this->present_index*this->dof_handler->get_fe().dofs_per_line+i]
    = index;
}





template <class DH>
void DoFObjectAccessor<1, DH>::set_vertex_dof_index (const unsigned int vertex,
                                                     const unsigned int i,
                                                     const unsigned int index,
					      const unsigned int fe_index) const

{
  Assert (fe_index == DoFHandler<1>::default_fe_index,
	  ExcMessage ("Only the default FE index is allowed for non-hp DoFHandler objects"));
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
  Assert (&this->dof_handler->get_fe() != 0,
	  typename BaseClass::ExcInvalidObject());
  Assert (vertex<2, ExcIndexRange (i,0,2));
  Assert (i<this->dof_handler->get_fe().dofs_per_vertex,
	  ExcIndexRange (i, 0, this->dof_handler->get_fe().dofs_per_vertex));

  const unsigned int dof_number = (this->vertex_index(vertex) *
				   this->dof_handler->get_fe().dofs_per_vertex +
				   i);
  this->dof_handler->vertex_dofs[dof_number] = index;
}



template <class DH>
template <class InputVector, typename number>
void
DoFObjectAccessor<1,DH>::get_dof_values (const InputVector &values,
                                         Vector<number>    &local_values) const
{
  typedef DoFAccessor<DH> BaseClass;

  Assert (dim==1, ExcInternalError());
  
  Assert (this->dof_handler != 0, typename BaseClass::ExcInvalidObject());
  Assert (&this->dof_handler->get_fe() != 0,
	  typename BaseClass::ExcInvalidObject());
  Assert (local_values.size() == this->dof_handler->get_fe().dofs_per_cell,
	  typename BaseClass::ExcVectorDoesNotMatch());
  Assert (values.size() == this->dof_handler->n_dofs(),
	  typename BaseClass::ExcVectorDoesNotMatch());
  Assert (this->has_children() == false,
	  typename BaseClass::ExcNotActive());
  
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
}



template <class DH>
template <class OutputVector, typename number>
void
DoFObjectAccessor<1,DH>::set_dof_values (const Vector<number> &local_values,
                                         OutputVector         &values) const
{
  typedef DoFAccessor<DH> BaseClass;

  Assert (dim==1, ExcInternalError());
  
  Assert (this->dof_handler != 0, typename BaseClass::ExcInvalidObject());
  Assert (&this->dof_handler->get_fe() != 0,
	  typename BaseClass::ExcInvalidObject());
  Assert (local_values.size() == this->dof_handler->get_fe().dofs_per_cell,
	  typename BaseClass::ExcVectorDoesNotMatch());
  Assert (values.size() == this->dof_handler->n_dofs(),
	  typename BaseClass::ExcVectorDoesNotMatch());
  Assert (this->has_children() == false,
	  typename BaseClass::ExcNotActive());
  
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
}



/*------------------------- Functions: DoFObjectAccessor<2,dim> -----------------------*/


template <class DH>
void DoFObjectAccessor<2, DH>::set_dof_index (const unsigned int i,
                                              const unsigned int index,
					      const unsigned int fe_index) const
{
  Assert (fe_index == DoFHandler<1>::default_fe_index,
	  ExcMessage ("Only the default FE index is allowed for non-hp DoFHandler objects"));
  typedef DoFAccessor<DH> BaseClass;

  Assert (this->dof_handler != 0,
	  typename BaseClass::ExcInvalidObject());
				   // make sure a FE has been selected
				   // and enough room was reserved
  Assert (&this->dof_handler->get_fe() != 0,
	  typename BaseClass::ExcInvalidObject());
  Assert (i<this->dof_handler->get_fe().dofs_per_quad,
	  ExcIndexRange (i, 0, this->dof_handler->get_fe().dofs_per_quad));

  this->dof_handler->levels[this->present_level]
    ->quad_dofs[this->present_index*this->dof_handler->get_fe().dofs_per_quad+i] = index;
}



template <class DH>
void
DoFObjectAccessor<2, DH>::set_vertex_dof_index (const unsigned int vertex,
                                                const unsigned int i,
                                                const unsigned int index,
					      const unsigned int fe_index) const
{
  Assert (fe_index == DoFHandler<1>::default_fe_index,
	  ExcMessage ("Only the default FE index is allowed for non-hp DoFHandler objects"));
  typedef DoFAccessor<DH> BaseClass;

  Assert (this->dof_handler != 0,
	  typename BaseClass::ExcInvalidObject());
  Assert (&this->dof_handler->get_fe() != 0,
	  typename BaseClass::ExcInvalidObject());
  Assert (vertex<4, ExcIndexRange (i,0,4));
  Assert (i<this->dof_handler->get_fe().dofs_per_vertex,
	  ExcIndexRange (i, 0, this->dof_handler->get_fe().dofs_per_vertex));

  const unsigned int dof_number = (this->vertex_index(vertex) *
				   this->dof_handler->get_fe().dofs_per_vertex +
				   i);
  this->dof_handler->vertex_dofs[dof_number] = index;
}



template <class DH>
template <class InputVector, typename number>
void
DoFObjectAccessor<2,DH>::get_dof_values (const InputVector &values,
                                         Vector<number>    &local_values) const
{
  typedef DoFAccessor<DH> BaseClass;

  Assert (dim==2, ExcInternalError());
  
  Assert (this->dof_handler != 0,
	  typename BaseClass::ExcInvalidObject());
  Assert (&this->dof_handler->get_fe() != 0,
	  typename BaseClass::ExcInvalidObject());
  Assert (local_values.size() == this->dof_handler->get_fe().dofs_per_cell,
	  typename BaseClass::ExcVectorDoesNotMatch());
  Assert (values.size() == this->dof_handler->n_dofs(),
	  typename BaseClass::ExcVectorDoesNotMatch());
  Assert (this->has_children() == false,
	  typename BaseClass::ExcNotActive());
  
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
}



template <class DH>
template <class OutputVector, typename number>
void
DoFObjectAccessor<2,DH>::set_dof_values (const Vector<number> &local_values,
                                         OutputVector         &values) const
{
  typedef DoFAccessor<DH> BaseClass;

  Assert (dim==2, ExcInternalError());
  
  Assert (this->dof_handler != 0,
	  typename BaseClass::ExcInvalidObject());
  Assert (&this->dof_handler->get_fe() != 0,
	  typename BaseClass::ExcInvalidObject());
  Assert (local_values.size() == this->dof_handler->get_fe().dofs_per_cell,
	  typename BaseClass::ExcVectorDoesNotMatch());
  Assert (values.size() == this->dof_handler->n_dofs(),
	  typename BaseClass::ExcVectorDoesNotMatch());
  Assert (this->has_children() == false,
	  typename BaseClass::ExcNotActive());
  
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
}


/*------------------------- Functions: DoFObjectAccessor<3,dim> -----------------------*/


template <class DH>
void DoFObjectAccessor<3, DH>::set_dof_index (const unsigned int i,
                                              const unsigned int index,
					      const unsigned int fe_index) const
{
  Assert (fe_index == DoFHandler<1>::default_fe_index,
	  ExcMessage ("Only the default FE index is allowed for non-hp DoFHandler objects"));
  typedef DoFAccessor<DH> BaseClass;

  Assert (this->dof_handler != 0,
	  typename BaseClass::ExcInvalidObject());
				   // make sure a FE has been selected
				   // and enough room was reserved
  Assert (&this->dof_handler->get_fe() != 0,
	  typename BaseClass::ExcInvalidObject());
  Assert (i<this->dof_handler->get_fe().dofs_per_hex,
	  ExcIndexRange (i, 0, this->dof_handler->get_fe().dofs_per_hex));

  this->dof_handler->levels[this->present_level]
    ->hex_dofs[this->present_index*this->dof_handler->get_fe().dofs_per_hex+i] = index;
}



template <class DH>
void DoFObjectAccessor<3, DH>::set_vertex_dof_index (const unsigned int vertex,
                                                     const unsigned int i,
                                                     const unsigned int index,
					      const unsigned int fe_index) const
{
  Assert (fe_index == DoFHandler<1>::default_fe_index,
	  ExcMessage ("Only the default FE index is allowed for non-hp DoFHandler objects"));
  typedef DoFAccessor<DH> BaseClass;

  Assert (this->dof_handler != 0,
	  typename BaseClass::ExcInvalidObject());
  Assert (&this->dof_handler->get_fe() != 0,
	  typename BaseClass::ExcInvalidObject());
  Assert (vertex<8,
	  ExcIndexRange (i,0,8));
  Assert (i<this->dof_handler->get_fe().dofs_per_vertex,
	  ExcIndexRange (i, 0, this->dof_handler->get_fe().dofs_per_vertex));

  const unsigned int dof_number = (this->vertex_index(vertex) *
				   this->dof_handler->get_fe().dofs_per_vertex +
				   i);
  this->dof_handler->vertex_dofs[dof_number] = index;
}



template <class DH>
template <class InputVector, typename number>
void
DoFObjectAccessor<3,DH>::get_dof_values (const InputVector &values,
                                         Vector<number>    &local_values) const
{
  typedef DoFAccessor<DH> BaseClass;

  Assert (dim==3, ExcInternalError());

  Assert (this->dof_handler != 0,
	  typename BaseClass::ExcInvalidObject());
  Assert (&this->dof_handler->get_fe() != 0, 
	  typename BaseClass::ExcInvalidObject());
  Assert (local_values.size() == this->dof_handler->get_fe().dofs_per_cell,
	  typename BaseClass::ExcVectorDoesNotMatch());
  Assert (values.size() == this->dof_handler->n_dofs(),
	  typename BaseClass::ExcVectorDoesNotMatch());
  Assert (this->has_children() == false,
	  typename BaseClass::ExcNotActive());
  
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
}



template <class DH>
template <class OutputVector, typename number>
void
DoFObjectAccessor<3,DH>::set_dof_values (const Vector<number> &local_values,
                                         OutputVector         &values) const
{
  typedef DoFAccessor<DH> BaseClass;

  Assert (dim==3, ExcInternalError());

  Assert (this->dof_handler != 0,
	  typename BaseClass::ExcInvalidObject());
  Assert (&this->dof_handler->get_fe() != 0,
	  typename BaseClass::ExcInvalidObject());
  Assert (local_values.size() == this->dof_handler->get_fe().dofs_per_cell,
	  typename BaseClass::ExcVectorDoesNotMatch());
  Assert (values.size() == this->dof_handler->n_dofs(),
	  typename BaseClass::ExcVectorDoesNotMatch());
  Assert (this->has_children() == false,
	  typename BaseClass::ExcNotActive());
  
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
}



// --------------- hp::DoFHandler specializations for 1d objects -----------


template <>
void DoFObjectAccessor<1, hp::DoFHandler<1> >::set_vertex_dof_index (const unsigned int /*vertex*/,
                                                                     const unsigned int /*i*/,
                                                                     const unsigned int /*index*/,
					      const unsigned int fe_index) const
{
  Assert (fe_index != hp::DoFHandler<1>::default_fe_index,
	  ExcMessage ("You need to specify a FE index when working with hp DoFHandlers"));
  Assert (false, ExcInternalError());
}

template <>
void DoFObjectAccessor<1, hp::DoFHandler<2> >::set_vertex_dof_index (const unsigned int /*vertex*/,
                                                                     const unsigned int /*i*/,
                                                                     const unsigned int /*index*/,
					      const unsigned int fe_index) const
{
  Assert (fe_index != hp::DoFHandler<2>::default_fe_index,
	  ExcMessage ("You need to specify a FE index when working with hp DoFHandlers"));
  Assert (false, ExcInternalError());
}

template <>
void DoFObjectAccessor<1, hp::DoFHandler<3> >::set_vertex_dof_index (const unsigned int /*vertex*/,
                                                                     const unsigned int /*i*/,
                                                                     const unsigned int /*index*/,
					      const unsigned int fe_index) const
{
  Assert (fe_index != hp::DoFHandler<3>::default_fe_index,
	  ExcMessage ("You need to specify a FE index when working with hp DoFHandlers"));
  Assert (false, ExcInternalError());
}


template <>
template <class InputVector, typename number>
void
DoFObjectAccessor<1,hp::DoFHandler<1> >::get_dof_values (const InputVector &/*values*/,
                                                         Vector<number>    &/*local_values*/) const
{
  Assert (false, ExcNotImplemented());
}

template <>
template <class InputVector, typename number>
void
DoFObjectAccessor<1,hp::DoFHandler<2> >::get_dof_values (const InputVector &/*values*/,
                                                         Vector<number>    &/*local_values*/) const
{
  Assert (false, ExcNotImplemented());
}

template <>
template <class InputVector, typename number>
void
DoFObjectAccessor<1,hp::DoFHandler<3> >::get_dof_values (const InputVector &/*values*/,
                                                         Vector<number>    &/*local_values*/) const
{
  Assert (false, ExcNotImplemented());
}


template <>
template <class OutputVector, typename number>
void
DoFObjectAccessor<1,hp::DoFHandler<1> >::set_dof_values (const Vector<number> &/*local_values*/,
                                                         OutputVector         &/*values*/) const
{
  Assert (false, ExcNotImplemented());
}

template <>
template <class OutputVector, typename number>
void
DoFObjectAccessor<1,hp::DoFHandler<2> >::set_dof_values (const Vector<number> &/*local_values*/,
                                                         OutputVector         &/*values*/) const
{
  Assert (false, ExcNotImplemented());
}

template <>
template <class OutputVector, typename number>
void
DoFObjectAccessor<1,hp::DoFHandler<3> >::set_dof_values (const Vector<number> &/*local_values*/,
                                                         OutputVector         &/*values*/) const
{
  Assert (false, ExcNotImplemented());
}




// --------------- hp::DoFHandler specializations for 2d objects -----------


template <>
void DoFObjectAccessor<2, hp::DoFHandler<2> >::set_vertex_dof_index (const unsigned int /*vertex*/,
                                                                     const unsigned int /*i*/,
                                                                     const unsigned int /*index*/,
					      const unsigned int fe_index) const
{
  Assert (fe_index != hp::DoFHandler<2>::default_fe_index,
	  ExcMessage ("You need to specify a FE index when working with hp DoFHandlers"));
  Assert (false, ExcInternalError());
}

template <>
void DoFObjectAccessor<2, hp::DoFHandler<3> >::set_vertex_dof_index (const unsigned int /*vertex*/,
                                                                     const unsigned int /*i*/,
                                                                     const unsigned int /*index*/,
					      const unsigned int fe_index) const
{
  Assert (fe_index != hp::DoFHandler<3>::default_fe_index,
	  ExcMessage ("You need to specify a FE index when working with hp DoFHandlers"));
  Assert (false, ExcInternalError());
}


template <>
template <class InputVector, typename number>
void
DoFObjectAccessor<2,hp::DoFHandler<2> >::get_dof_values (const InputVector &/*values*/,
                                                         Vector<number>    &/*local_values*/) const
{
  Assert (false, ExcNotImplemented());
}

template <>
template <class InputVector, typename number>
void
DoFObjectAccessor<2,hp::DoFHandler<3> >::get_dof_values (const InputVector &/*values*/,
                                                         Vector<number>    &/*local_values*/) const
{
  Assert (false, ExcNotImplemented());
}


template <>
template <class OutputVector, typename number>
void
DoFObjectAccessor<2,hp::DoFHandler<2> >::set_dof_values (const Vector<number> &/*local_values*/,
                                                         OutputVector         &/*values*/) const
{
  Assert (false, ExcNotImplemented());
}

template <>
template <class OutputVector, typename number>
void
DoFObjectAccessor<2,hp::DoFHandler<3> >::set_dof_values (const Vector<number> &/*local_values*/,
                                                         OutputVector         &/*values*/) const
{
  Assert (false, ExcNotImplemented());
}




// --------------- hp::DoFHandler specializations for 3d objects -----------


template <>
void DoFObjectAccessor<3, hp::DoFHandler<3> >::set_vertex_dof_index (const unsigned int /*vertex*/,
                                                                     const unsigned int /*i*/,
                                                                     const unsigned int /*index*/,
					      const unsigned int fe_index) const
{
  Assert (fe_index != hp::DoFHandler<3>::default_fe_index,
	  ExcMessage ("You need to specify a FE index when working with hp DoFHandlers"));
  Assert (false, ExcInternalError());
}

template <>
template <class InputVector, typename number>
void
DoFObjectAccessor<3,hp::DoFHandler<3> >::get_dof_values (const InputVector &/*values*/,
                                                         Vector<number>    &/*local_values*/) const
{
  Assert (false, ExcNotImplemented());
}

template <>
template <class OutputVector, typename number>
void
DoFObjectAccessor<3,hp::DoFHandler<3> >::set_dof_values (const Vector<number> &/*local_values*/,
                                                         OutputVector         &/*values*/) const
{
  Assert (false, ExcNotImplemented());
}



/*------------------------- Functions: DoFCellAccessor -----------------------*/


#if deal_II_dimension == 1

template <>
TriaIterator<1, DoFObjectAccessor<0,DoFHandler<1> > >
DoFCellAccessor<DoFHandler<1> >::face (const unsigned int) const
{
  Assert (false, ExcImpossibleInDim(1));
  return TriaIterator<1, DoFObjectAccessor<0,DoFHandler<1> > >();
}


template <>
void
DoFCellAccessor<DoFHandler<1> >::update_cell_dof_indices_cache () const
{
  Assert (static_cast<unsigned int>(this->present_level) < this->dof_handler->levels.size(),
          ExcMessage ("DoFHandler not initialized"));
  
  Assert (this->dof_handler != 0, DoFAccessor<DoFHandler<1> >::ExcInvalidObject());
  Assert (&this->get_fe() != 0, DoFAccessor<DoFHandler<1> >::ExcInvalidObject());

				   // check as in documentation that
				   // cell is either active, or dofs
				   // are only in vertices. otherwise
				   // simply don't update the cache at
				   // all. the get_dof_indices
				   // function will then make sure we
				   // don't access the invalid data
  if (this->has_children()
      &&
      (this->get_fe().dofs_per_cell !=
       this->get_fe().dofs_per_vertex * GeometryInfo<1>::vertices_per_cell))
    return;
  
  const unsigned int dofs_per_vertex = this->get_fe().dofs_per_vertex,
		     dofs_per_line   = this->get_fe().dofs_per_line,
		     dofs_per_cell   = this->get_fe().dofs_per_cell;

				   // make sure the cache is at least
				   // as big as we need it when
				   // writing to the last element of
				   // this cell
  Assert (this->present_index * dofs_per_cell + dofs_per_cell
	  <=
	  this->dof_handler->levels[this->present_level]
	  ->cell_dof_indices_cache.size(),
	  ExcInternalError());

  std::vector<unsigned int>::iterator next
    = (this->dof_handler->levels[this->present_level]
       ->cell_dof_indices_cache.begin() + this->present_index * dofs_per_cell);
  
  for (unsigned int vertex=0; vertex<2; ++vertex)
    for (unsigned int d=0; d<dofs_per_vertex; ++d)
      *next++ = vertex_dof_index(vertex,d);
  for (unsigned int d=0; d<dofs_per_line; ++d)
    *next++ = dof_index(d);
}



template <>
void
DoFCellAccessor<hp::DoFHandler<1> >::update_cell_dof_indices_cache () const
{
				   // not implemented, but should also
				   // not be called
  Assert (false, ExcNotImplemented());
}



#endif


#if deal_II_dimension == 2

template <>
TriaIterator<2, DoFObjectAccessor<1,DoFHandler<2> > >
DoFCellAccessor<DoFHandler<2> >::face (const unsigned int i) const
{
  return this->line(i);
}



template <>
void
DoFCellAccessor<DoFHandler<2> >::update_cell_dof_indices_cache () const
{
  Assert (static_cast<unsigned int>(this->present_level) < this->dof_handler->levels.size(),
          ExcMessage ("DoFHandler not initialized"));

  Assert (this->dof_handler != 0, DoFAccessor<DoFHandler<2> >::ExcInvalidObject());
  Assert (&this->get_fe() != 0, DoFAccessor<DoFHandler<2> >::ExcInvalidObject());

				   // check as in documentation that
				   // cell is either active, or dofs
				   // are only in vertices. otherwise
				   // simply don't update the cache at
				   // all. the get_dof_indices
				   // function will then make sure we
				   // don't access the invalid data
  if (this->has_children()
      &&
      (this->get_fe().dofs_per_cell !=
       this->get_fe().dofs_per_vertex * GeometryInfo<2>::vertices_per_cell))
    return;
  
  const unsigned int dofs_per_vertex = this->get_fe().dofs_per_vertex,
		     dofs_per_line   = this->get_fe().dofs_per_line,
		     dofs_per_quad   = this->get_fe().dofs_per_quad,
		     dofs_per_cell   = this->get_fe().dofs_per_cell;

				   // make sure the cache is at least
				   // as big as we need it when
				   // writing to the last element of
				   // this cell
  Assert (this->present_index * dofs_per_cell + dofs_per_cell
	  <=
	  this->dof_handler->levels[this->present_level]
	  ->cell_dof_indices_cache.size(),
	  ExcInternalError());

  std::vector<unsigned int>::iterator next
    = (this->dof_handler->levels[this->present_level]
       ->cell_dof_indices_cache.begin() + this->present_index * dofs_per_cell);

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
void
DoFCellAccessor<hp::DoFHandler<2> >::update_cell_dof_indices_cache () const
{
				   // not implemented, but should also
				   // not be called
  Assert (false, ExcNotImplemented());
}


#endif


#if deal_II_dimension == 3

template <>
TriaIterator<3, DoFObjectAccessor<2,DoFHandler<3> > >
DoFCellAccessor<DoFHandler<3> >::face (const unsigned int i) const
{
  return this->quad(i);
}



template <>
void
DoFCellAccessor<DoFHandler<3> >::update_cell_dof_indices_cache () const
{
  Assert (static_cast<unsigned int>(this->present_level) < this->dof_handler->levels.size(),
          ExcMessage ("DoFHandler not initialized"));

  Assert (this->dof_handler != 0, DoFAccessor<DoFHandler<3> >::ExcInvalidObject());
  Assert (&this->get_fe() != 0, DoFAccessor<DoFHandler<3> >::ExcInvalidObject());

				   // check as in documentation that
				   // cell is either active, or dofs
				   // are only in vertices. otherwise
				   // simply don't update the cache at
				   // all. the get_dof_indices
				   // function will then make sure we
				   // don't access the invalid data
  if (this->has_children()
      &&
      (this->get_fe().dofs_per_cell !=
       this->get_fe().dofs_per_vertex * GeometryInfo<3>::vertices_per_cell))
    return;
  
  const unsigned int dofs_per_vertex = this->get_fe().dofs_per_vertex,
		     dofs_per_line   = this->get_fe().dofs_per_line,
		     dofs_per_quad   = this->get_fe().dofs_per_quad,
		     dofs_per_hex    = this->get_fe().dofs_per_hex,
		     dofs_per_cell   = this->get_fe().dofs_per_cell;

				   // make sure the cache is at least
				   // as big as we need it when
				   // writing to the last element of
				   // this cell
  Assert (this->present_index * dofs_per_cell + dofs_per_cell
	  <=
	  this->dof_handler->levels[this->present_level]
	  ->cell_dof_indices_cache.size(),
	  ExcInternalError());

  std::vector<unsigned int>::iterator next
    = (this->dof_handler->levels[this->present_level]
       ->cell_dof_indices_cache.begin() + this->present_index * dofs_per_cell);

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


template <>
void
DoFCellAccessor<hp::DoFHandler<3> >::update_cell_dof_indices_cache () const
{
				   // not implemented, but should also
				   // not be called
  Assert (false, ExcNotImplemented());
}


#endif


template <class DH>
TriaIterator<DH::dimension,DoFCellAccessor<DH> >
DoFCellAccessor<DH>::neighbor_child_on_subface (const unsigned int face,
						const unsigned int subface) const
{
  const TriaIterator<dim,CellAccessor<dim> > q
    = CellAccessor<dim>::neighbor_child_on_subface (face, subface);
  return TriaIterator<dim,DoFCellAccessor<DH> > (this->tria,
						 q->level (),
						 q->index (),
						 this->dof_handler);
}



template <class DH>
template <class InputVector, typename number>
void
DoFCellAccessor<DH>::
get_interpolated_dof_values (const InputVector &values,
			     Vector<number>    &interpolated_values) const
{
  const FiniteElement<dim> &fe            = this->get_fe();
  const unsigned int        dofs_per_cell = fe.dofs_per_cell;
  
  typedef DoFAccessor<DH> BaseClass;
  
  Assert (this->dof_handler != 0,
	  typename BaseClass::ExcInvalidObject());
  Assert (&fe != 0,
	  typename BaseClass::ExcInvalidObject());
  Assert (interpolated_values.size() == dofs_per_cell,
	  typename BaseClass::ExcVectorDoesNotMatch());
  Assert (values.size() == this->dof_handler->n_dofs(),
	  typename BaseClass::ExcVectorDoesNotMatch());

  if (!this->has_children())
				     // if this cell has no children: simply
				     // return the exact values on this cell
    this->get_dof_values (values, interpolated_values);
  else
				     // otherwise clobber them from the children
    {
      Vector<number> tmp1(dofs_per_cell);
      Vector<number> tmp2(dofs_per_cell);
      
      interpolated_values = 0;

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
                                       // an example for the latter
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
                                       // individually. in fact, there
                                       // are even weird finite
                                       // elements (for example the
                                       // Raviart-Thomas element)
                                       // which have shape functions
                                       // that are additive (interior
                                       // ones) and others that are
                                       // overwriting (face degrees of
                                       // freedom that need to be
                                       // continuous across the
                                       // face). to avoid checking
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
	  fe.get_restriction_matrix(child).vmult (tmp2, tmp1);

                                           // and add up or set them
                                           // in the output vector
	  for (unsigned int i=0; i<dofs_per_cell; ++i)
            if (restriction_is_additive[i]) 
              interpolated_values(i) += tmp2(i);
            else
              if (tmp2(i) != 0)
                interpolated_values(i) = tmp2(i);
	}
    }
}



template <class DH>
template <class OutputVector, typename number>
void
DoFCellAccessor<DH>::
set_dof_values_by_interpolation (const Vector<number> &local_values,
				 OutputVector         &values) const
{
  const unsigned int dofs_per_cell = this->get_fe().dofs_per_cell;
  typedef DoFAccessor<DH> BaseClass;
  
  Assert (this->dof_handler != 0,
	  typename BaseClass::ExcInvalidObject());
  Assert (&this->get_fe() != 0,
	  typename BaseClass::ExcInvalidObject());
  Assert (local_values.size() == dofs_per_cell,
	  typename BaseClass::ExcVectorDoesNotMatch());
  Assert (values.size() == this->dof_handler->n_dofs(),
	  typename BaseClass::ExcVectorDoesNotMatch());

  if (!this->has_children())
                                     // if this cell has no children: simply
				     // set the values on this cell
    this->set_dof_values (local_values, values);
  else
				     // otherwise distribute them to the children
    {
      Vector<number> tmp(dofs_per_cell);
      
      for (unsigned int child=0; child<GeometryInfo<dim>::children_per_cell;
	   ++child)
	{
					   // prolong the given data
					   // to the present cell
	  this->get_fe().get_prolongation_matrix(child)
            .vmult (tmp, local_values);
	  this->child(child)->set_dof_values_by_interpolation (tmp, values);
	}
    }
}



// --------------------------------------------------------------------------
// explicit instantiations (for DoFHandler)


template
void
DoFObjectAccessor<1,DoFHandler<deal_II_dimension> >::get_dof_values<Vector<double>,double>
(const Vector<double>&, Vector<double>&) const;
template
void
DoFObjectAccessor<1,DoFHandler<deal_II_dimension> >::set_dof_values<Vector<double>,double>
(const Vector<double>&, Vector<double>&) const;
template
void
DoFObjectAccessor<1,DoFHandler<deal_II_dimension> >::get_dof_values<Vector<double>,float>
(const Vector<double>&, Vector<float>&) const;
template
void
DoFObjectAccessor<1,DoFHandler<deal_II_dimension> >::set_dof_values<Vector<float>,double>
(const Vector<double>&, Vector<float>&) const;
template
void
DoFObjectAccessor<1,DoFHandler<deal_II_dimension> >::get_dof_values<Vector<float>,double>
(const Vector<float>&, Vector<double>&) const;
template
void
DoFObjectAccessor<1,DoFHandler<deal_II_dimension> >::set_dof_values<Vector<double>,float>
(const Vector<float>&, Vector<double>&) const;
template
void
DoFObjectAccessor<1,DoFHandler<deal_II_dimension> >::get_dof_values<Vector<float>,float>
(const Vector<float>&, Vector<float>&) const;
template
void
DoFObjectAccessor<1,DoFHandler<deal_II_dimension> >::set_dof_values<Vector<float>,float>
(const Vector<float>&, Vector<float>&) const;


// for block vector
template
void
DoFObjectAccessor<1,DoFHandler<deal_II_dimension> >::get_dof_values<BlockVector<double>,double>
(const BlockVector<double> &, Vector<double>&) const;
template
void
DoFObjectAccessor<1,DoFHandler<deal_II_dimension> >::set_dof_values<BlockVector<double>,double>
(const Vector<double>&, BlockVector<double>&) const;
template
void
DoFObjectAccessor<1,DoFHandler<deal_II_dimension> >::get_dof_values<BlockVector<double>,float>
(const BlockVector<double> &, Vector<float>&) const;
template
void
DoFObjectAccessor<1,DoFHandler<deal_II_dimension> >::set_dof_values<BlockVector<float>,double>
(const Vector<double>&, BlockVector<float>&) const;
template
void
DoFObjectAccessor<1,DoFHandler<deal_II_dimension> >::get_dof_values<BlockVector<float>,double>
(const BlockVector<float> &, Vector<double>&) const;
template
void
DoFObjectAccessor<1,DoFHandler<deal_II_dimension> >::set_dof_values<BlockVector<double>,float>
(const Vector<float>&, BlockVector<double>&) const;
template
void
DoFObjectAccessor<1,DoFHandler<deal_II_dimension> >::get_dof_values<BlockVector<float>,float>
(const BlockVector<float>&, Vector<float>&) const;
template
void
DoFObjectAccessor<1,DoFHandler<deal_II_dimension> >::set_dof_values<BlockVector<float>,float>
(const Vector<float>&, BlockVector<float>&) const;

// for Petsc vectors
#ifdef DEAL_II_USE_PETSC
template
void
DoFObjectAccessor<1,DoFHandler<deal_II_dimension> >::get_dof_values<PETScWrappers::Vector,double>
(const PETScWrappers::Vector &, Vector<double>&) const;
template
void
DoFObjectAccessor<1,DoFHandler<deal_II_dimension> >::get_dof_values<PETScWrappers::Vector,float>
(const PETScWrappers::Vector &, Vector<float>&) const;
template
void
DoFObjectAccessor<1,DoFHandler<deal_II_dimension> >::set_dof_values<PETScWrappers::Vector,double>
(const Vector<double> &, PETScWrappers::Vector&) const;
template
void
DoFObjectAccessor<1,DoFHandler<deal_II_dimension> >::set_dof_values<PETScWrappers::Vector,float>
(const Vector<float>&, PETScWrappers::Vector&) const;

template
void
DoFObjectAccessor<1,DoFHandler<deal_II_dimension> >::get_dof_values<PETScWrappers::BlockVector,double>
(const PETScWrappers::BlockVector &, Vector<double>&) const;
template
void
DoFObjectAccessor<1,DoFHandler<deal_II_dimension> >::get_dof_values<PETScWrappers::BlockVector,float>
(const PETScWrappers::BlockVector &, Vector<float>&) const;
template
void
DoFObjectAccessor<1,DoFHandler<deal_II_dimension> >::set_dof_values<PETScWrappers::BlockVector,double>
(const Vector<double> &, PETScWrappers::BlockVector&) const;
template
void
DoFObjectAccessor<1,DoFHandler<deal_II_dimension> >::set_dof_values<PETScWrappers::BlockVector,float>
(const Vector<float>&, PETScWrappers::BlockVector&) const;
#endif

#if deal_II_dimension >= 2
template
void
DoFObjectAccessor<2,DoFHandler<deal_II_dimension> >::get_dof_values<Vector<double>,double>
(const Vector<double>&, Vector<double>&) const;
template
void
DoFObjectAccessor<2,DoFHandler<deal_II_dimension> >::set_dof_values<Vector<double>,double>
(const Vector<double>&, Vector<double>&) const;
template
void
DoFObjectAccessor<2,DoFHandler<deal_II_dimension> >::get_dof_values<Vector<double>,float>
(const Vector<double>&, Vector<float>&) const;
template
void
DoFObjectAccessor<2,DoFHandler<deal_II_dimension> >::set_dof_values<Vector<float>,double>
(const Vector<double>&, Vector<float>&) const;
template
void
DoFObjectAccessor<2,DoFHandler<deal_II_dimension> >::get_dof_values<Vector<float>,double>
(const Vector<float>&, Vector<double>&) const;
template
void
DoFObjectAccessor<2,DoFHandler<deal_II_dimension> >::set_dof_values<Vector<double>,float>
(const Vector<float>&, Vector<double>&) const;
template
void
DoFObjectAccessor<2,DoFHandler<deal_II_dimension> >::get_dof_values<Vector<float>,float>
(const Vector<float>&, Vector<float>&) const;
template
void
DoFObjectAccessor<2,DoFHandler<deal_II_dimension> >::set_dof_values<Vector<float>,float>
(const Vector<float>&, Vector<float>&) const;


// for block vector
template
void
DoFObjectAccessor<2,DoFHandler<deal_II_dimension> >::get_dof_values<BlockVector<double>,double>
(const BlockVector<double>&, Vector<double>&) const;
template
void
DoFObjectAccessor<2,DoFHandler<deal_II_dimension> >::set_dof_values<BlockVector<double>,double>
(const Vector<double>&, BlockVector<double>&) const;
template
void
DoFObjectAccessor<2,DoFHandler<deal_II_dimension> >::get_dof_values<BlockVector<double>,float>
(const BlockVector<double>&, Vector<float>&) const;
template
void
DoFObjectAccessor<2,DoFHandler<deal_II_dimension> >::set_dof_values<BlockVector<float>,double>
(const Vector<double>&, BlockVector<float>&) const;
template
void
DoFObjectAccessor<2,DoFHandler<deal_II_dimension> >::get_dof_values<BlockVector<float>,double>
(const BlockVector<float>&, Vector<double>&) const;
template
void
DoFObjectAccessor<2,DoFHandler<deal_II_dimension> >::set_dof_values<BlockVector<double>,float>
(const Vector<float>&, BlockVector<double>&) const;
template
void
DoFObjectAccessor<2,DoFHandler<deal_II_dimension> >::get_dof_values<BlockVector<float>,float>
(const BlockVector<float>&, Vector<float>&) const;
template
void
DoFObjectAccessor<2,DoFHandler<deal_II_dimension> >::set_dof_values<BlockVector<float>,float>
(const Vector<float>&, BlockVector<float>&) const;

// for Petsc vectors
#ifdef DEAL_II_USE_PETSC
template
void
DoFObjectAccessor<2,DoFHandler<deal_II_dimension> >::get_dof_values<PETScWrappers::Vector,double>
(const PETScWrappers::Vector &, Vector<double>&) const;
template
void
DoFObjectAccessor<2,DoFHandler<deal_II_dimension> >::get_dof_values<PETScWrappers::Vector,float>
(const PETScWrappers::Vector &, Vector<float>&) const;
template
void
DoFObjectAccessor<2,DoFHandler<deal_II_dimension> >::set_dof_values<PETScWrappers::Vector,double>
(const Vector<double> &, PETScWrappers::Vector&) const;
template
void
DoFObjectAccessor<2,DoFHandler<deal_II_dimension> >::set_dof_values<PETScWrappers::Vector,float>
(const Vector<float>&, PETScWrappers::Vector&) const;

template
void
DoFObjectAccessor<2,DoFHandler<deal_II_dimension> >::get_dof_values<PETScWrappers::BlockVector,double>
(const PETScWrappers::BlockVector &, Vector<double>&) const;
template
void
DoFObjectAccessor<2,DoFHandler<deal_II_dimension> >::get_dof_values<PETScWrappers::BlockVector,float>
(const PETScWrappers::BlockVector &, Vector<float>&) const;
template
void
DoFObjectAccessor<2,DoFHandler<deal_II_dimension> >::set_dof_values<PETScWrappers::BlockVector,double>
(const Vector<double> &, PETScWrappers::BlockVector&) const;
template
void
DoFObjectAccessor<2,DoFHandler<deal_II_dimension> >::set_dof_values<PETScWrappers::BlockVector,float>
(const Vector<float>&, PETScWrappers::BlockVector&) const;
#endif
#endif

#if deal_II_dimension >= 3
template
void
DoFObjectAccessor<3,DoFHandler<deal_II_dimension> >::get_dof_values<Vector<double>,double>
(const Vector<double>&, Vector<double>&) const;
template
void
DoFObjectAccessor<3,DoFHandler<deal_II_dimension> >::set_dof_values<Vector<double>,double>
(const Vector<double>&, Vector<double>&) const;
template
void
DoFObjectAccessor<3,DoFHandler<deal_II_dimension> >::get_dof_values<Vector<double>,float>
(const Vector<double>&, Vector<float>&) const;
template
void
DoFObjectAccessor<3,DoFHandler<deal_II_dimension> >::set_dof_values<Vector<float>,double>
(const Vector<double>&, Vector<float>&) const;
template
void
DoFObjectAccessor<3,DoFHandler<deal_II_dimension> >::get_dof_values<Vector<float>,double>
(const Vector<float>&, Vector<double>&) const;
template
void
DoFObjectAccessor<3,DoFHandler<deal_II_dimension> >::set_dof_values<Vector<double>,float>
(const Vector<float>&, Vector<double>&) const;
template
void
DoFObjectAccessor<3,DoFHandler<deal_II_dimension> >::get_dof_values<Vector<float>,float>
(const Vector<float>&, Vector<float>&) const;
template
void
DoFObjectAccessor<3,DoFHandler<deal_II_dimension> >::set_dof_values<Vector<float>,float>
(const Vector<float>&, Vector<float>&) const;



// for block vector
template
void
DoFObjectAccessor<3,DoFHandler<deal_II_dimension> >::get_dof_values<BlockVector<double>,double>
(const BlockVector<double>&, Vector<double>&) const;
template
void
DoFObjectAccessor<3,DoFHandler<deal_II_dimension> >::set_dof_values<BlockVector<double>,double>
(const Vector<double>&, BlockVector<double>&) const;
template
void
DoFObjectAccessor<3,DoFHandler<deal_II_dimension> >::get_dof_values<BlockVector<double>,float>
(const BlockVector<double>&, Vector<float>&) const;
template
void
DoFObjectAccessor<3,DoFHandler<deal_II_dimension> >::set_dof_values<BlockVector<float>,double>
(const Vector<double>&, BlockVector<float>&) const;
template
void
DoFObjectAccessor<3,DoFHandler<deal_II_dimension> >::get_dof_values<BlockVector<float>,double>
(const BlockVector<float>&, Vector<double>&) const;
template
void
DoFObjectAccessor<3,DoFHandler<deal_II_dimension> >::set_dof_values<BlockVector<double>,float>
(const Vector<float>&, BlockVector<double>&) const;
template
void
DoFObjectAccessor<3,DoFHandler<deal_II_dimension> >::get_dof_values<BlockVector<float>,float>
(const BlockVector<float>&, Vector<float>&) const;
template
void
DoFObjectAccessor<3,DoFHandler<deal_II_dimension> >::set_dof_values<BlockVector<float>,float>
(const Vector<float>&, BlockVector<float>&) const;

// for Petsc vectors
#if DEAL_II_USE_PETSC
template
void
DoFObjectAccessor<3,DoFHandler<deal_II_dimension> >::get_dof_values<PETScWrappers::Vector,double>
(const PETScWrappers::Vector &, Vector<double>&) const;
template
void
DoFObjectAccessor<3,DoFHandler<deal_II_dimension> >::get_dof_values<PETScWrappers::Vector,float>
(const PETScWrappers::Vector &, Vector<float>&) const;
template
void
DoFObjectAccessor<3,DoFHandler<deal_II_dimension> >::set_dof_values<PETScWrappers::Vector,double>
(const Vector<double> &, PETScWrappers::Vector&) const;
template
void
DoFObjectAccessor<3,DoFHandler<deal_II_dimension> >::set_dof_values<PETScWrappers::Vector,float>
(const Vector<float>&, PETScWrappers::Vector&) const;

template
void
DoFObjectAccessor<3,DoFHandler<deal_II_dimension> >::get_dof_values<PETScWrappers::BlockVector,double>
(const PETScWrappers::BlockVector &, Vector<double>&) const;
template
void
DoFObjectAccessor<3,DoFHandler<deal_II_dimension> >::get_dof_values<PETScWrappers::BlockVector,float>
(const PETScWrappers::BlockVector &, Vector<float>&) const;
template
void
DoFObjectAccessor<3,DoFHandler<deal_II_dimension> >::set_dof_values<PETScWrappers::BlockVector,double>
(const Vector<double> &, PETScWrappers::BlockVector&) const;
template
void
DoFObjectAccessor<3,DoFHandler<deal_II_dimension> >::set_dof_values<PETScWrappers::BlockVector,float>
(const Vector<float>&, PETScWrappers::BlockVector&) const;
#endif
#endif



template
void
DoFCellAccessor<DoFHandler<deal_II_dimension> >::
get_interpolated_dof_values<Vector<double>,double>
(const Vector<double>&, Vector<double>&) const;
template
void
DoFCellAccessor<DoFHandler<deal_II_dimension> >::
set_dof_values_by_interpolation<Vector<double>,double>
(const Vector<double>&, Vector<double>&) const;

template
void
DoFCellAccessor<DoFHandler<deal_II_dimension> >::
get_interpolated_dof_values<Vector<double>,float>
(const Vector<double>&, Vector<float>&) const;
template
void
DoFCellAccessor<DoFHandler<deal_II_dimension> >::
set_dof_values_by_interpolation<Vector<float>,double>
(const Vector<double>&, Vector<float>&) const;

template
void
DoFCellAccessor<DoFHandler<deal_II_dimension> >::
get_interpolated_dof_values<Vector<float>,double>
(const Vector<float>&, Vector<double>&) const;
template
void
DoFCellAccessor<DoFHandler<deal_II_dimension> >::
set_dof_values_by_interpolation<Vector<double>,float>
(const Vector<float>&, Vector<double>&) const;

template
void
DoFCellAccessor<DoFHandler<deal_II_dimension> >::
get_interpolated_dof_values<Vector<float>,float>
(const Vector<float>&, Vector<float>&) const;
template
void
DoFCellAccessor<DoFHandler<deal_II_dimension> >::
set_dof_values_by_interpolation<Vector<float>,float>
(const Vector<float>&, Vector<float>&) const;


template
void
DoFCellAccessor<DoFHandler<deal_II_dimension> >::
get_interpolated_dof_values<BlockVector<double>,double>
(const BlockVector<double>&, Vector<double>&) const;

template
void
DoFCellAccessor<DoFHandler<deal_II_dimension> >::
set_dof_values_by_interpolation<BlockVector<double>,double>
(const Vector<double>&, BlockVector<double>&) const;

template
void
DoFCellAccessor<DoFHandler<deal_II_dimension> >::
get_interpolated_dof_values<BlockVector<double>,float>
(const BlockVector<double>&, Vector<float>&) const;

template
void
DoFCellAccessor<DoFHandler<deal_II_dimension> >::
set_dof_values_by_interpolation<BlockVector<float>,double>
(const Vector<double>&, BlockVector<float>&) const;

template
void
DoFCellAccessor<DoFHandler<deal_II_dimension> >::
get_interpolated_dof_values<BlockVector<float>,double>
(const BlockVector<float>&, Vector<double>&) const;
template
void
DoFCellAccessor<DoFHandler<deal_II_dimension> >::
set_dof_values_by_interpolation<BlockVector<double>,float>
(const Vector<float>&, BlockVector<double>&) const;

template
void
DoFCellAccessor<DoFHandler<deal_II_dimension> >::
get_interpolated_dof_values<BlockVector<float>,float>
(const BlockVector<float>&, Vector<float>&) const;
template
void
DoFCellAccessor<DoFHandler<deal_II_dimension> >::
set_dof_values_by_interpolation<BlockVector<float>,float>
(const Vector<float>&, BlockVector<float>&) const;

// for Petsc vectors
#if DEAL_II_USE_PETSC
template
void
DoFCellAccessor<DoFHandler<deal_II_dimension> >::
get_interpolated_dof_values<PETScWrappers::Vector,double>
(const PETScWrappers::Vector&, Vector<double>&) const;
template
void
DoFCellAccessor<DoFHandler<deal_II_dimension> >::
set_dof_values_by_interpolation<PETScWrappers::Vector,double>
(const Vector<double>&, PETScWrappers::Vector&) const;

template
void
DoFCellAccessor<DoFHandler<deal_II_dimension> >::
get_interpolated_dof_values<PETScWrappers::Vector,float>
(const PETScWrappers::Vector&, Vector<float>&) const;
template
void
DoFCellAccessor<DoFHandler<deal_II_dimension> >::
set_dof_values_by_interpolation<PETScWrappers::Vector,float>
(const Vector<float>&, PETScWrappers::Vector&) const;


template
void
DoFCellAccessor<DoFHandler<deal_II_dimension> >::
get_interpolated_dof_values<PETScWrappers::BlockVector,double>
(const PETScWrappers::BlockVector&, Vector<double>&) const;
template
void
DoFCellAccessor<DoFHandler<deal_II_dimension> >::
set_dof_values_by_interpolation<PETScWrappers::BlockVector,double>
(const Vector<double>&, PETScWrappers::BlockVector&) const;

template
void
DoFCellAccessor<DoFHandler<deal_II_dimension> >::
get_interpolated_dof_values<PETScWrappers::BlockVector,float>
(const PETScWrappers::BlockVector&, Vector<float>&) const;
template
void
DoFCellAccessor<DoFHandler<deal_II_dimension> >::
set_dof_values_by_interpolation<PETScWrappers::BlockVector,float>
(const Vector<float>&, PETScWrappers::BlockVector&) const;

#endif


template class DoFAccessor<DoFHandler<deal_II_dimension> >;

#if deal_II_dimension == 1
template class DoFObjectAccessor<1, DoFHandler<1> >;
#endif

#if deal_II_dimension == 2
template class DoFObjectAccessor<1, DoFHandler<2> >;
template class DoFObjectAccessor<2, DoFHandler<2> >;

template class TriaRawIterator   <2,DoFObjectAccessor<1, DoFHandler<2> > >;
template class TriaIterator      <2,DoFObjectAccessor<1, DoFHandler<2> > >;
template class TriaActiveIterator<2,DoFObjectAccessor<1, DoFHandler<2> > >;
#endif

#if deal_II_dimension == 3
template class DoFObjectAccessor<1, DoFHandler<3> >;
template class DoFObjectAccessor<2, DoFHandler<3> >;
template class DoFObjectAccessor<3, DoFHandler<3> >;

template class TriaRawIterator   <3,DoFObjectAccessor<1, DoFHandler<3> > >;
template class TriaIterator      <3,DoFObjectAccessor<1, DoFHandler<3> > >;
template class TriaActiveIterator<3,DoFObjectAccessor<1, DoFHandler<3> > >;
template class TriaRawIterator   <3,DoFObjectAccessor<2, DoFHandler<3> > >;
template class TriaIterator      <3,DoFObjectAccessor<2, DoFHandler<3> > >;
template class TriaActiveIterator<3,DoFObjectAccessor<2, DoFHandler<3> > >;
#endif


template class DoFCellAccessor<DoFHandler<deal_II_dimension> >;

template class TriaRawIterator   <deal_II_dimension,DoFCellAccessor<DoFHandler<deal_II_dimension> > >;
template class TriaIterator      <deal_II_dimension,DoFCellAccessor<DoFHandler<deal_II_dimension> > >;
template class TriaActiveIterator<deal_II_dimension,DoFCellAccessor<DoFHandler<deal_II_dimension> > >;


// --------------------------------------------------------------------------
// explicit instantiations (for hp::DoFHandler)


template
void
DoFObjectAccessor<1,hp::DoFHandler<deal_II_dimension> >::get_dof_values<Vector<double>,double>
(const Vector<double>&, Vector<double>&) const;
template
void
DoFObjectAccessor<1,hp::DoFHandler<deal_II_dimension> >::set_dof_values<Vector<double>,double>
(const Vector<double>&, Vector<double>&) const;
template
void
DoFObjectAccessor<1,hp::DoFHandler<deal_II_dimension> >::get_dof_values<Vector<double>,float>
(const Vector<double>&, Vector<float>&) const;
template
void
DoFObjectAccessor<1,hp::DoFHandler<deal_II_dimension> >::set_dof_values<Vector<float>,double>
(const Vector<double>&, Vector<float>&) const;
template
void
DoFObjectAccessor<1,hp::DoFHandler<deal_II_dimension> >::get_dof_values<Vector<float>,double>
(const Vector<float>&, Vector<double>&) const;
template
void
DoFObjectAccessor<1,hp::DoFHandler<deal_II_dimension> >::set_dof_values<Vector<double>,float>
(const Vector<float>&, Vector<double>&) const;
template
void
DoFObjectAccessor<1,hp::DoFHandler<deal_II_dimension> >::get_dof_values<Vector<float>,float>
(const Vector<float>&, Vector<float>&) const;
template
void
DoFObjectAccessor<1,hp::DoFHandler<deal_II_dimension> >::set_dof_values<Vector<float>,float>
(const Vector<float>&, Vector<float>&) const;


// for block vector
template
void
DoFObjectAccessor<1,hp::DoFHandler<deal_II_dimension> >::get_dof_values<BlockVector<double>,double>
(const BlockVector<double> &, Vector<double>&) const;
template
void
DoFObjectAccessor<1,hp::DoFHandler<deal_II_dimension> >::set_dof_values<BlockVector<double>,double>
(const Vector<double>&, BlockVector<double>&) const;
template
void
DoFObjectAccessor<1,hp::DoFHandler<deal_II_dimension> >::get_dof_values<BlockVector<double>,float>
(const BlockVector<double> &, Vector<float>&) const;
template
void
DoFObjectAccessor<1,hp::DoFHandler<deal_II_dimension> >::set_dof_values<BlockVector<float>,double>
(const Vector<double>&, BlockVector<float>&) const;
template
void
DoFObjectAccessor<1,hp::DoFHandler<deal_II_dimension> >::get_dof_values<BlockVector<float>,double>
(const BlockVector<float> &, Vector<double>&) const;
template
void
DoFObjectAccessor<1,hp::DoFHandler<deal_II_dimension> >::set_dof_values<BlockVector<double>,float>
(const Vector<float>&, BlockVector<double>&) const;
template
void
DoFObjectAccessor<1,hp::DoFHandler<deal_II_dimension> >::get_dof_values<BlockVector<float>,float>
(const BlockVector<float>&, Vector<float>&) const;
template
void
DoFObjectAccessor<1,hp::DoFHandler<deal_II_dimension> >::set_dof_values<BlockVector<float>,float>
(const Vector<float>&, BlockVector<float>&) const;

// for Petsc vectors
#if DEAL_II_USE_PETSC
template
void
DoFObjectAccessor<1,hp::DoFHandler<deal_II_dimension> >::get_dof_values<PETScWrappers::Vector,double>
(const PETScWrappers::Vector &, Vector<double>&) const;
template
void
DoFObjectAccessor<1,hp::DoFHandler<deal_II_dimension> >::get_dof_values<PETScWrappers::Vector,float>
(const PETScWrappers::Vector &, Vector<float>&) const;
template
void
DoFObjectAccessor<1,hp::DoFHandler<deal_II_dimension> >::set_dof_values<PETScWrappers::Vector,double>
(const Vector<double> &, PETScWrappers::Vector&) const;
template
void
DoFObjectAccessor<1,hp::DoFHandler<deal_II_dimension> >::set_dof_values<PETScWrappers::Vector,float>
(const Vector<float>&, PETScWrappers::Vector&) const;
#endif


#if deal_II_dimension >= 2
template
void
DoFObjectAccessor<2,hp::DoFHandler<deal_II_dimension> >::get_dof_values<Vector<double>,double>
(const Vector<double>&, Vector<double>&) const;
template
void
DoFObjectAccessor<2,hp::DoFHandler<deal_II_dimension> >::set_dof_values<Vector<double>,double>
(const Vector<double>&, Vector<double>&) const;
template
void
DoFObjectAccessor<2,hp::DoFHandler<deal_II_dimension> >::get_dof_values<Vector<double>,float>
(const Vector<double>&, Vector<float>&) const;
template
void
DoFObjectAccessor<2,hp::DoFHandler<deal_II_dimension> >::set_dof_values<Vector<float>,double>
(const Vector<double>&, Vector<float>&) const;
template
void
DoFObjectAccessor<2,hp::DoFHandler<deal_II_dimension> >::get_dof_values<Vector<float>,double>
(const Vector<float>&, Vector<double>&) const;
template
void
DoFObjectAccessor<2,hp::DoFHandler<deal_II_dimension> >::set_dof_values<Vector<double>,float>
(const Vector<float>&, Vector<double>&) const;
template
void
DoFObjectAccessor<2,hp::DoFHandler<deal_II_dimension> >::get_dof_values<Vector<float>,float>
(const Vector<float>&, Vector<float>&) const;
template
void
DoFObjectAccessor<2,hp::DoFHandler<deal_II_dimension> >::set_dof_values<Vector<float>,float>
(const Vector<float>&, Vector<float>&) const;


// for block vector
template
void
DoFObjectAccessor<2,hp::DoFHandler<deal_II_dimension> >::get_dof_values<BlockVector<double>,double>
(const BlockVector<double>&, Vector<double>&) const;
template
void
DoFObjectAccessor<2,hp::DoFHandler<deal_II_dimension> >::set_dof_values<BlockVector<double>,double>
(const Vector<double>&, BlockVector<double>&) const;
template
void
DoFObjectAccessor<2,hp::DoFHandler<deal_II_dimension> >::get_dof_values<BlockVector<double>,float>
(const BlockVector<double>&, Vector<float>&) const;
template
void
DoFObjectAccessor<2,hp::DoFHandler<deal_II_dimension> >::set_dof_values<BlockVector<float>,double>
(const Vector<double>&, BlockVector<float>&) const;
template
void
DoFObjectAccessor<2,hp::DoFHandler<deal_II_dimension> >::get_dof_values<BlockVector<float>,double>
(const BlockVector<float>&, Vector<double>&) const;
template
void
DoFObjectAccessor<2,hp::DoFHandler<deal_II_dimension> >::set_dof_values<BlockVector<double>,float>
(const Vector<float>&, BlockVector<double>&) const;
template
void
DoFObjectAccessor<2,hp::DoFHandler<deal_II_dimension> >::get_dof_values<BlockVector<float>,float>
(const BlockVector<float>&, Vector<float>&) const;
template
void
DoFObjectAccessor<2,hp::DoFHandler<deal_II_dimension> >::set_dof_values<BlockVector<float>,float>
(const Vector<float>&, BlockVector<float>&) const;

// for Petsc vectors
#if DEAL_II_USE_PETSC
template
void
DoFObjectAccessor<2,hp::DoFHandler<deal_II_dimension> >::get_dof_values<PETScWrappers::Vector,double>
(const PETScWrappers::Vector &, Vector<double>&) const;
template
void
DoFObjectAccessor<2,hp::DoFHandler<deal_II_dimension> >::get_dof_values<PETScWrappers::Vector,float>
(const PETScWrappers::Vector &, Vector<float>&) const;
template
void
DoFObjectAccessor<2,hp::DoFHandler<deal_II_dimension> >::set_dof_values<PETScWrappers::Vector,double>
(const Vector<double> &, PETScWrappers::Vector&) const;
template
void
DoFObjectAccessor<2,hp::DoFHandler<deal_II_dimension> >::set_dof_values<PETScWrappers::Vector,float>
(const Vector<float>&, PETScWrappers::Vector&) const;
#endif

#endif

#if deal_II_dimension >= 3
// for double
template
void
DoFObjectAccessor<3,hp::DoFHandler<deal_II_dimension> >::get_dof_values<Vector<double>,double>
(const Vector<double>&, Vector<double>&) const;
template
void
DoFObjectAccessor<3,hp::DoFHandler<deal_II_dimension> >::set_dof_values<Vector<double>,double>
(const Vector<double>&, Vector<double>&) const;
template
void
DoFObjectAccessor<3,hp::DoFHandler<deal_II_dimension> >::get_dof_values<Vector<double>,float>
(const Vector<double>&, Vector<float>&) const;
template
void
DoFObjectAccessor<3,hp::DoFHandler<deal_II_dimension> >::set_dof_values<Vector<float>,double>
(const Vector<double>&, Vector<float>&) const;
template
void
DoFObjectAccessor<3,hp::DoFHandler<deal_II_dimension> >::get_dof_values<Vector<float>,double>
(const Vector<float>&, Vector<double>&) const;
template
void
DoFObjectAccessor<3,hp::DoFHandler<deal_II_dimension> >::set_dof_values<Vector<double>,float>
(const Vector<float>&, Vector<double>&) const;
template
void
DoFObjectAccessor<3,hp::DoFHandler<deal_II_dimension> >::get_dof_values<Vector<float>,float>
(const Vector<float>&, Vector<float>&) const;
template
void
DoFObjectAccessor<3,hp::DoFHandler<deal_II_dimension> >::set_dof_values<Vector<float>,float>
(const Vector<float>&, Vector<float>&) const;



// for block vector
template
void
DoFObjectAccessor<3,hp::DoFHandler<deal_II_dimension> >::get_dof_values<BlockVector<double>,double>
(const BlockVector<double>&, Vector<double>&) const;
template
void
DoFObjectAccessor<3,hp::DoFHandler<deal_II_dimension> >::set_dof_values<BlockVector<double>,double>
(const Vector<double>&, BlockVector<double>&) const;
template
void
DoFObjectAccessor<3,hp::DoFHandler<deal_II_dimension> >::get_dof_values<BlockVector<double>,float>
(const BlockVector<double>&, Vector<float>&) const;
template
void
DoFObjectAccessor<3,hp::DoFHandler<deal_II_dimension> >::set_dof_values<BlockVector<float>,double>
(const Vector<double>&, BlockVector<float>&) const;
template
void
DoFObjectAccessor<3,hp::DoFHandler<deal_II_dimension> >::get_dof_values<BlockVector<float>,double>
(const BlockVector<float>&, Vector<double>&) const;
template
void
DoFObjectAccessor<3,hp::DoFHandler<deal_II_dimension> >::set_dof_values<BlockVector<double>,float>
(const Vector<float>&, BlockVector<double>&) const;
template
void
DoFObjectAccessor<3,hp::DoFHandler<deal_II_dimension> >::get_dof_values<BlockVector<float>,float>
(const BlockVector<float>&, Vector<float>&) const;
template
void
DoFObjectAccessor<3,hp::DoFHandler<deal_II_dimension> >::set_dof_values<BlockVector<float>,float>
(const Vector<float>&, BlockVector<float>&) const;


// for Petsc vectors
#ifdef DEAL_II_USE_PETSC
template
void
DoFObjectAccessor<3,hp::DoFHandler<deal_II_dimension> >::get_dof_values<PETScWrappers::Vector,double>
(const PETScWrappers::Vector &, Vector<double>&) const;
template
void
DoFObjectAccessor<3,hp::DoFHandler<deal_II_dimension> >::get_dof_values<PETScWrappers::Vector,float>
(const PETScWrappers::Vector &, Vector<float>&) const;
template
void
DoFObjectAccessor<3,hp::DoFHandler<deal_II_dimension> >::set_dof_values<PETScWrappers::Vector,double>
(const Vector<double> &, PETScWrappers::Vector&) const;
template
void
DoFObjectAccessor<3,hp::DoFHandler<deal_II_dimension> >::set_dof_values<PETScWrappers::Vector,float>
(const Vector<float>&, PETScWrappers::Vector&) const;

template
void
DoFObjectAccessor<3,hp::DoFHandler<deal_II_dimension> >::get_dof_values<PETScWrappers::BlockVector,double>
(const PETScWrappers::BlockVector &, Vector<double>&) const;
template
void
DoFObjectAccessor<3,hp::DoFHandler<deal_II_dimension> >::get_dof_values<PETScWrappers::BlockVector,float>
(const PETScWrappers::BlockVector &, Vector<float>&) const;
template
void
DoFObjectAccessor<3,hp::DoFHandler<deal_II_dimension> >::set_dof_values<PETScWrappers::BlockVector,double>
(const Vector<double>&, PETScWrappers::BlockVector &) const;
template
void
DoFObjectAccessor<3,hp::DoFHandler<deal_II_dimension> >::set_dof_values<PETScWrappers::BlockVector,float>
(const Vector<float>&, PETScWrappers::BlockVector&) const;
#endif
#endif



template
void
DoFCellAccessor<hp::DoFHandler<deal_II_dimension> >::
get_interpolated_dof_values<Vector<double>,double>
(const Vector<double>&, Vector<double>&) const;
template
void
DoFCellAccessor<hp::DoFHandler<deal_II_dimension> >::
set_dof_values_by_interpolation<Vector<double>,double>
(const Vector<double>&, Vector<double>&) const;

template
void
DoFCellAccessor<hp::DoFHandler<deal_II_dimension> >::
get_interpolated_dof_values<Vector<double>,float>
(const Vector<double>&, Vector<float>&) const;
template
void
DoFCellAccessor<hp::DoFHandler<deal_II_dimension> >::
set_dof_values_by_interpolation<Vector<float>,double>
(const Vector<double>&, Vector<float>&) const;

template
void
DoFCellAccessor<hp::DoFHandler<deal_II_dimension> >::
get_interpolated_dof_values<Vector<float>,double>
(const Vector<float>&, Vector<double>&) const;
template
void
DoFCellAccessor<hp::DoFHandler<deal_II_dimension> >::
set_dof_values_by_interpolation<Vector<double>,float>
(const Vector<float>&, Vector<double>&) const;

template
void
DoFCellAccessor<hp::DoFHandler<deal_II_dimension> >::
get_interpolated_dof_values<Vector<float>,float>
(const Vector<float>&, Vector<float>&) const;
template
void
DoFCellAccessor<hp::DoFHandler<deal_II_dimension> >::
set_dof_values_by_interpolation<Vector<float>,float>
(const Vector<float>&, Vector<float>&) const;


template
void
DoFCellAccessor<hp::DoFHandler<deal_II_dimension> >::
get_interpolated_dof_values<BlockVector<double>,double>
(const BlockVector<double>&, Vector<double>&) const;

template
void
DoFCellAccessor<hp::DoFHandler<deal_II_dimension> >::
set_dof_values_by_interpolation<BlockVector<double>,double>
(const Vector<double>&, BlockVector<double>&) const;

template
void
DoFCellAccessor<hp::DoFHandler<deal_II_dimension> >::
get_interpolated_dof_values<BlockVector<double>,float>
(const BlockVector<double>&, Vector<float>&) const;

template
void
DoFCellAccessor<hp::DoFHandler<deal_II_dimension> >::
set_dof_values_by_interpolation<BlockVector<float>,double>
(const Vector<double>&, BlockVector<float>&) const;

template
void
DoFCellAccessor<hp::DoFHandler<deal_II_dimension> >::
get_interpolated_dof_values<BlockVector<float>,double>
(const BlockVector<float>&, Vector<double>&) const;
template
void
DoFCellAccessor<hp::DoFHandler<deal_II_dimension> >::
set_dof_values_by_interpolation<BlockVector<double>,float>
(const Vector<float>&, BlockVector<double>&) const;

template
void
DoFCellAccessor<hp::DoFHandler<deal_II_dimension> >::
get_interpolated_dof_values<BlockVector<float>,float>
(const BlockVector<float>&, Vector<float>&) const;
template
void
DoFCellAccessor<hp::DoFHandler<deal_II_dimension> >::
set_dof_values_by_interpolation<BlockVector<float>,float>
(const Vector<float>&, BlockVector<float>&) const;

// for Petsc vectors
#ifdef DEAL_II_USE_PETSC
template
void
DoFCellAccessor<hp::DoFHandler<deal_II_dimension> >::
get_interpolated_dof_values<PETScWrappers::Vector,double>
(const PETScWrappers::Vector&, Vector<double>&) const;
template
void
DoFCellAccessor<hp::DoFHandler<deal_II_dimension> >::
set_dof_values_by_interpolation<PETScWrappers::Vector,double>
(const Vector<double>&, PETScWrappers::Vector&) const;

template
void
DoFCellAccessor<hp::DoFHandler<deal_II_dimension> >::
get_interpolated_dof_values<PETScWrappers::Vector,float>
(const PETScWrappers::Vector&, Vector<float>&) const;
template
void
DoFCellAccessor<hp::DoFHandler<deal_II_dimension> >::
set_dof_values_by_interpolation<PETScWrappers::Vector,float>
(const Vector<float>&, PETScWrappers::Vector&) const;


template
void
DoFCellAccessor<hp::DoFHandler<deal_II_dimension> >::
get_interpolated_dof_values<PETScWrappers::BlockVector,double>
(const PETScWrappers::BlockVector&, Vector<double>&) const;
template
void
DoFCellAccessor<hp::DoFHandler<deal_II_dimension> >::
set_dof_values_by_interpolation<PETScWrappers::BlockVector,double>
(const Vector<double>&, PETScWrappers::BlockVector&) const;

template
void
DoFCellAccessor<hp::DoFHandler<deal_II_dimension> >::
get_interpolated_dof_values<PETScWrappers::BlockVector,float>
(const PETScWrappers::BlockVector&, Vector<float>&) const;
template
void
DoFCellAccessor<hp::DoFHandler<deal_II_dimension> >::
set_dof_values_by_interpolation<PETScWrappers::BlockVector,float>
(const Vector<float>&, PETScWrappers::BlockVector&) const;
#endif


template class DoFAccessor<hp::DoFHandler<deal_II_dimension> >;

#if deal_II_dimension == 1
template class DoFObjectAccessor<1, hp::DoFHandler<1> >;
#endif

#if deal_II_dimension == 2
template class DoFObjectAccessor<1, hp::DoFHandler<2> >;
template class DoFObjectAccessor<2, hp::DoFHandler<2> >;

template class TriaRawIterator   <2,DoFObjectAccessor<1, hp::DoFHandler<2> > >;
template class TriaIterator      <2,DoFObjectAccessor<1, hp::DoFHandler<2> > >;
template class TriaActiveIterator<2,DoFObjectAccessor<1, hp::DoFHandler<2> > >;
#endif


#if deal_II_dimension == 3
template class DoFObjectAccessor<1, hp::DoFHandler<3> >;
template class DoFObjectAccessor<2, hp::DoFHandler<3> >;
template class DoFObjectAccessor<3, hp::DoFHandler<3> >;

template class TriaRawIterator   <3,DoFObjectAccessor<1, hp::DoFHandler<3> > >;
template class TriaIterator      <3,DoFObjectAccessor<1, hp::DoFHandler<3> > >;
template class TriaActiveIterator<3,DoFObjectAccessor<1, hp::DoFHandler<3> > >;
template class TriaRawIterator   <3,DoFObjectAccessor<2, hp::DoFHandler<3> > >;
template class TriaIterator      <3,DoFObjectAccessor<2, hp::DoFHandler<3> > >;
template class TriaActiveIterator<3,DoFObjectAccessor<2, hp::DoFHandler<3> > >;
#endif


template class DoFCellAccessor<hp::DoFHandler<deal_II_dimension> >;

template class TriaRawIterator   <deal_II_dimension,DoFCellAccessor<hp::DoFHandler<deal_II_dimension> > >;
template class TriaIterator      <deal_II_dimension,DoFCellAccessor<hp::DoFHandler<deal_II_dimension> > >;
template class TriaActiveIterator<deal_II_dimension,DoFCellAccessor<hp::DoFHandler<deal_II_dimension> > >;

