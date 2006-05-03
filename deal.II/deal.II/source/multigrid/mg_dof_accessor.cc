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


#include <grid/tria_iterator_base.h>
#include <dofs/dof_levels.h>
#include <multigrid/mg_dof_accessor.h>
#include <multigrid/mg_dof_handler.h>
#include <grid/tria_iterator.h>
#include <grid/tria_iterator.templates.h>
#include <fe/fe.h>

#include <lac/vector.h>
#include <lac/full_matrix.h>
#include <lac/sparse_matrix.h>

//TODO:[GK] Inline simple functions in 1d and 3d

/* ------------------------ MGDoFAccessor --------------------------- */


template <int structdim, int dim>
MGDoFAccessor<structdim, dim>::MGDoFAccessor ()
		:
		BaseClass (0,
			   deal_II_numbers::invalid_unsigned_int,
			   deal_II_numbers::invalid_unsigned_int,
			   0),
		mg_dof_handler(0)
{
  Assert (false, ExcInvalidObject());
}



template <int structdim, int dim>
MGDoFAccessor<structdim, dim>::MGDoFAccessor (const Triangulation<dim> *tria,
					      const int                 level,
					      const int                 index,
					      const AccessorData       *local_data)
		:
		BaseClass (tria, level, index, local_data),
		mg_dof_handler (const_cast<MGDoFHandler<dim>*>(local_data))
{}



template <int structdim, int dim>
void
MGDoFAccessor<structdim, dim>::set_mg_dof_handler (MGDoFHandler<dim> *dh)
{
  typedef DoFAccessor<dim,DoFHandler<dim> > BaseClass;
  Assert (dh != 0, typename BaseClass::ExcInvalidObject());
  mg_dof_handler = dh;
}



template <int structdim, int dim>
MGDoFAccessor<structdim, dim> &
MGDoFAccessor<structdim, dim>::operator = (const MGDoFAccessor &da)
{
  BaseClass::operator= (*this);
  
  set_dof_handler (da.mg_dof_handler);
  return *this;
}



template <int structdim,int dim>
unsigned int
MGDoFAccessor<structdim, dim>::mg_vertex_dof_index (const unsigned int vertex,
						    const unsigned int i) const
{
  Assert (this->dof_handler != 0,
	  typename BaseClass::ExcInvalidObject());
  Assert (this->mg_dof_handler != 0,
	  typename BaseClass::ExcInvalidObject());
  Assert (&this->dof_handler->get_fe() != 0,
	  typename BaseClass::ExcInvalidObject());
  Assert (vertex<GeometryInfo<structdim>::vertices_per_cell,
	  ExcIndexRange (i,0,GeometryInfo<structdim>::vertices_per_cell));
  Assert (i<this->dof_handler->get_fe().dofs_per_vertex,
	  ExcIndexRange (i, 0, this->dof_handler->get_fe().dofs_per_vertex));

  return (this->mg_dof_handler->mg_vertex_dofs[this->vertex_index(vertex)]
	  .get_index (this->present_level, i, this->dof_handler->get_fe().dofs_per_vertex));
}


template <int structdim, int dim>
void
MGDoFAccessor<structdim, dim>::set_mg_vertex_dof_index (const unsigned int vertex,
							const unsigned int i,
							const unsigned int index) const
{
  Assert (this->dof_handler != 0,
	  typename BaseClass::ExcInvalidObject());
  Assert (this->mg_dof_handler != 0,
	  typename BaseClass::ExcInvalidObject());
  Assert (&this->dof_handler->get_fe() != 0,
	  typename BaseClass::ExcInvalidObject());
  Assert (vertex<GeometryInfo<structdim>::vertices_per_cell,
	  ExcIndexRange (i,0,GeometryInfo<structdim>::vertices_per_cell));
  Assert (i<this->dof_handler->get_fe().dofs_per_vertex,
	  ExcIndexRange (i, 0, this->dof_handler->get_fe().dofs_per_vertex));

  this->mg_dof_handler->mg_vertex_dofs[this->vertex_index(vertex)]
    .set_index (this->present_level, i, this->dof_handler->get_fe().dofs_per_vertex, index);
}



template <int structdim, int dim>
unsigned int
MGDoFAccessor<structdim, dim>::mg_dof_index (const unsigned int i) const
{
  return this->mg_dof_handler->mg_levels[this->present_level]
    ->get_dof_index(*this->dof_handler,
		    this->present_index,
		    0,
		    i,
		    internal::StructuralDimension<structdim>());
}


template <int structdim, int dim>
void MGDoFAccessor<structdim, dim>::set_mg_dof_index (const unsigned int i,
						    const unsigned int index) const
{
  this->mg_dof_handler->mg_levels[this->present_level]
    ->set_dof_index(*this->dof_handler,
		    this->present_index,
		    0,
		    i,
		    index,
		    internal::StructuralDimension<structdim>());
}



template <int structdim, int dim>
TriaIterator<dim,MGDoFObjectAccessor<structdim, dim> >
MGDoFAccessor<structdim, dim>::child (const unsigned int i) const
{
  TriaIterator<dim,MGDoFObjectAccessor<structdim, dim> > q (this->tria,
							    this->present_level+1,
							    this->child_index (i),
							    this->mg_dof_handler);
  
				   // make sure that we either created
				   // a past-the-end iterator or one
				   // pointing to a used cell
  Assert ((q.state() == IteratorState::past_the_end)
	  ||
	  q->used(),
	  typename TriaAccessor<dim>::ExcUnusedCellAsChild());

  return q;
}



template <int structdim, int dim>
void
MGDoFAccessor<structdim, dim>::copy_from (const MGDoFAccessor &a)
{
  DoFObjectAccessor<structdim, DoFHandler<dim> >::copy_from (a);
  this->set_mg_dof_handler (a.mg_dof_handler);
}



/* ------------------------ MGDoFObjectAccessor<0> --------------------------- */

template <int dim>
MGDoFObjectAccessor<0,dim>::MGDoFObjectAccessor (const Triangulation<dim> *,
						 const int                 ,
						 const int                 ,
						 const AccessorData       *)
{
  Assert (false, ExcInternalError());
}
  


/* ------------------------ MGDoFObjectAccessor<0> --------------------------- */


template <int dim>
MGDoFObjectAccessor<1, dim>::MGDoFObjectAccessor (const Triangulation<dim> *tria,
						  const int                 level,
						  const int                 index,
						  const AccessorData       *local_data)
		:
		MGDoFAccessor<1,dim> (tria,level,index,local_data)
{}



template <int dim>
void
MGDoFObjectAccessor<1, dim>::get_mg_dof_indices (std::vector<unsigned int> &dof_indices) const
{
  Assert (this->dof_handler != 0,
	  typename BaseClass::ExcInvalidObject());
  Assert (this->mg_dof_handler != 0,
	  typename BaseClass::ExcInvalidObject());
  Assert (&this->dof_handler->get_fe() != 0,
	  typename BaseClass::ExcInvalidObject());
  Assert (dof_indices.size() == (2*this->dof_handler->get_fe().dofs_per_vertex +
				 this->dof_handler->get_fe().dofs_per_line),
	  typename BaseClass::ExcVectorDoesNotMatch());

  const unsigned int dofs_per_vertex = this->dof_handler->get_fe().dofs_per_vertex,
		     dofs_per_line   = this->dof_handler->get_fe().dofs_per_line;
  std::vector<unsigned int>::iterator next = dof_indices.begin();
  for (unsigned int vertex=0; vertex<2; ++vertex)
    for (unsigned int d=0; d<dofs_per_vertex; ++d)
      *next++ = this->mg_vertex_dof_index(vertex,d);
  for (unsigned int d=0; d<dofs_per_line; ++d)
    *next++ = this->mg_dof_index(d);
  
  Assert (next == dof_indices.end(),
	  ExcInternalError());
}


template <int dim>
template <typename number>
void
MGDoFObjectAccessor<1,dim>::get_mg_dof_values (const Vector<number> &values,
					       Vector<number>       &dof_values) const
{
  Assert (this->dof_handler != 0,
	  typename BaseClass::ExcInvalidObject());
  Assert (this->mg_dof_handler != 0,
	  typename BaseClass::ExcInvalidObject());
  Assert (&this->dof_handler->get_fe() != 0,
	  typename BaseClass::ExcInvalidObject());
  Assert (dof_values.size() == this->dof_handler->get_fe().dofs_per_cell,
	  typename BaseClass::ExcVectorDoesNotMatch());
  Assert (values.size() == this->dof_handler->n_dofs(),
	  typename BaseClass::ExcVectorDoesNotMatch());

  const unsigned int dofs_per_vertex = this->dof_handler->get_fe().dofs_per_vertex,
		     dofs_per_line   = this->dof_handler->get_fe().dofs_per_line;
  typename Vector<number>::iterator next_dof_value=dof_values.begin();
  for (unsigned int vertex=0; vertex<2; ++vertex)
    for (unsigned int d=0; d<dofs_per_vertex; ++d)
      *next_dof_value++ = values(this->mg_vertex_dof_index(vertex,d));
  for (unsigned int d=0; d<dofs_per_line; ++d)
    *next_dof_value++ = values(this->mg_dof_index(d));
  
  Assert (next_dof_value == dof_values.end(),
	  ExcInternalError());
}



/* ------------------------ MGDoFObjectAccessor<2> --------------------------- */

template <int dim>
MGDoFObjectAccessor<2, dim>::MGDoFObjectAccessor (const Triangulation<dim> *tria,
						  const int                 level,
						  const int                 index,
						  const AccessorData       *local_data)
		:
		MGDoFAccessor<2,dim> (tria,level,index,local_data)
{}




template <int dim>
void
MGDoFObjectAccessor<2, dim>::get_mg_dof_indices (std::vector<unsigned int> &dof_indices) const
{
  Assert (this->dof_handler != 0,
	  typename BaseClass::ExcInvalidObject());
  Assert (this->mg_dof_handler != 0,
	  typename BaseClass::ExcInvalidObject());
  Assert (&this->dof_handler->get_fe() != 0,
	  typename BaseClass::ExcInvalidObject());
  Assert (dof_indices.size() == (4*this->dof_handler->get_fe().dofs_per_vertex +
				 4*this->dof_handler->get_fe().dofs_per_line +
				 this->dof_handler->get_fe().dofs_per_quad),
	  typename BaseClass::ExcVectorDoesNotMatch());

  const unsigned int dofs_per_vertex = this->dof_handler->get_fe().dofs_per_vertex,
		     dofs_per_line   = this->dof_handler->get_fe().dofs_per_line,
		     dofs_per_quad   = this->dof_handler->get_fe().dofs_per_quad;
  std::vector<unsigned int>::iterator next = dof_indices.begin();
  for (unsigned int vertex=0; vertex<4; ++vertex)
    for (unsigned int d=0; d<dofs_per_vertex; ++d)
      *next++ = this->mg_vertex_dof_index(vertex,d);
  for (unsigned int line=0; line<4; ++line)
    for (unsigned int d=0; d<dofs_per_line; ++d)
      *next++ = this->line(line)->mg_dof_index(d);
  for (unsigned int d=0; d<dofs_per_quad; ++d)
    *next++ = this->mg_dof_index(d);
  
  Assert (next == dof_indices.end(),
	  ExcInternalError());
}


template <int dim>
template <typename number>
void
MGDoFObjectAccessor<2,dim>::get_mg_dof_values (const Vector<number> &values,
					       Vector<number>       &dof_values) const
{
  Assert (this->dof_handler != 0,
	  typename BaseClass::ExcInvalidObject());
  Assert (this->mg_dof_handler != 0,
	  typename BaseClass::ExcInvalidObject());
  Assert (&this->dof_handler->get_fe() != 0,
	  typename BaseClass::ExcInvalidObject());
  Assert (dof_values.size() == this->dof_handler->get_fe().dofs_per_cell,
	  typename BaseClass::ExcVectorDoesNotMatch());
  Assert (values.size() == this->mg_dof_handler->n_dofs(this->present_level),
	  typename BaseClass::ExcVectorDoesNotMatch());

  const unsigned int dofs_per_vertex = this->dof_handler->get_fe().dofs_per_vertex,
		     dofs_per_line   = this->dof_handler->get_fe().dofs_per_line,
		     dofs_per_quad   = this->dof_handler->get_fe().dofs_per_quad;
  typename Vector<number>::iterator next_dof_value=dof_values.begin();
  for (unsigned int vertex=0; vertex<4; ++vertex)
    for (unsigned int d=0; d<dofs_per_vertex; ++d)
      *next_dof_value++ = values(this->mg_vertex_dof_index(vertex,d));
  for (unsigned int line=0; line<4; ++line)
    for (unsigned int d=0; d<dofs_per_line; ++d)
      *next_dof_value++ = values(this->line(line)->mg_dof_index(d));
  for (unsigned int d=0; d<dofs_per_quad; ++d)
    *next_dof_value++ = values(this->mg_dof_index(d));
  
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


/* ------------------------ MGDoFObjectAccessor<3> --------------------------- */

template <int dim>
MGDoFObjectAccessor<3, dim>::MGDoFObjectAccessor (const Triangulation<dim> *tria,
						  const int                 level,
						  const int                 index,
						  const AccessorData       *local_data)
		:
		MGDoFAccessor<3,dim> (tria,level,index,local_data)
{}




template <int dim>
void
MGDoFObjectAccessor<3, dim>::get_mg_dof_indices (std::vector<unsigned int> &dof_indices) const
{
  Assert (this->dof_handler != 0,
	  typename BaseClass::ExcInvalidObject());
  Assert (this->mg_dof_handler != 0,
	  typename BaseClass::ExcInvalidObject());
  Assert (&this->dof_handler->get_fe() != 0,
	  typename BaseClass::ExcInvalidObject());
  Assert (dof_indices.size() == (8*this->dof_handler->get_fe().dofs_per_vertex +
				 12*this->dof_handler->get_fe().dofs_per_line +
				 6*this->dof_handler->get_fe().dofs_per_quad +
				 this->dof_handler->get_fe().dofs_per_hex),
	  typename BaseClass::ExcVectorDoesNotMatch());

  const unsigned int dofs_per_vertex = this->dof_handler->get_fe().dofs_per_vertex,
		     dofs_per_line   = this->dof_handler->get_fe().dofs_per_line,
		     dofs_per_quad   = this->dof_handler->get_fe().dofs_per_quad,
		     dofs_per_hex    = this->dof_handler->get_fe().dofs_per_hex;
  std::vector<unsigned int>::iterator next = dof_indices.begin();
  for (unsigned int vertex=0; vertex<8; ++vertex)
    for (unsigned int d=0; d<dofs_per_vertex; ++d)
      *next++ = this->mg_vertex_dof_index(vertex,d);
  for (unsigned int line=0; line<12; ++line)
    for (unsigned int d=0; d<dofs_per_line; ++d)
      *next++ = this->line(line)->mg_dof_index(d);
  for (unsigned int quad=0; quad<6; ++quad)
    for (unsigned int d=0; d<dofs_per_quad; ++d)
      *next++ = this->quad(quad)->mg_dof_index(d);
  for (unsigned int d=0; d<dofs_per_hex; ++d)
    *next++ = this->mg_dof_index(d);
  
  Assert (next == dof_indices.end(),
	  ExcInternalError());
}


template <int dim>
template <typename number>
void
MGDoFObjectAccessor<3,dim>::get_mg_dof_values (const Vector<number> &values,
					       Vector<number>       &dof_values) const
{
  Assert (this->dof_handler != 0,
	  typename BaseClass::ExcInvalidObject());
  Assert (this->mg_dof_handler != 0,
	  typename BaseClass::ExcInvalidObject());
  Assert (&this->dof_handler->get_fe() != 0,
	  typename BaseClass::ExcInvalidObject());
  Assert (dof_values.size() == this->dof_handler->get_fe().dofs_per_cell,
	  typename BaseClass::ExcVectorDoesNotMatch());
  Assert (values.size() == this->mg_dof_handler->n_dofs(this->present_level),
	  typename BaseClass::ExcVectorDoesNotMatch());

  const unsigned int dofs_per_vertex = this->dof_handler->get_fe().dofs_per_vertex,
		     dofs_per_line   = this->dof_handler->get_fe().dofs_per_line,
		     dofs_per_quad   = this->dof_handler->get_fe().dofs_per_quad,
		     dofs_per_hex    = this->dof_handler->get_fe().dofs_per_hex;
  typename Vector<number>::iterator next_dof_value=dof_values.begin();
  for (unsigned int vertex=0; vertex<8; ++vertex)
    for (unsigned int d=0; d<dofs_per_vertex; ++d)
      *next_dof_value++ = values(this->mg_vertex_dof_index(vertex,d));
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
MGDoFObjectAccessor<3, dim>::line (const unsigned int i) const
{
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
MGDoFObjectAccessor<3, dim>::quad (const unsigned int i) const
{
  Assert (i<12, ExcIndexRange (i, 0, 6));

  return TriaIterator<dim,MGDoFObjectAccessor<2, dim> >
    (
      this->tria,
      this->present_level,
      this->quad_index (i),
      this->mg_dof_handler
    );
}



/*------------------------- Functions: MGDoFCellAccessor -----------------------*/

template <int dim>
MGDoFCellAccessor<dim>::MGDoFCellAccessor (const Triangulation<dim> *tria,
					   const int                 level,
					   const int                 index,
					   const AccessorData       *local_data)
		:
		MGDoFObjectAccessor<dim, dim> (tria,level,index,local_data)
{}


template <int dim>
TriaIterator<dim,MGDoFCellAccessor<dim> >
MGDoFCellAccessor<dim>::neighbor (const unsigned int i) const
{
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
MGDoFCellAccessor<dim>::child (const unsigned int i) const
{
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
  Assert (false, ExcImpossibleInDim(1));
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



template <int dim>
TriaIterator<dim,MGDoFCellAccessor<dim> >
MGDoFCellAccessor<dim>::
neighbor_child_on_subface (const unsigned int face,
                           const unsigned int subface) const
{
  const TriaIterator<dim,CellAccessor<dim> > q
    = CellAccessor<dim>::neighbor_child_on_subface (face, subface);
  return TriaIterator<dim,MGDoFCellAccessor<dim> > (this->tria,
						    q->level (),
						    q->index (),
						    this->mg_dof_handler);
}




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
template class MGDoFAccessor<1, 1>;

template class MGDoFObjectAccessor<1, 1>;
template class MGDoFCellAccessor<1>;

template class TriaRawIterator<1,MGDoFCellAccessor<1> >;
template class TriaIterator<1,MGDoFCellAccessor<1> >;
template class TriaActiveIterator<1,MGDoFCellAccessor<1> >;
#endif

#if deal_II_dimension == 2
template class MGDoFAccessor<1, 2>;
template class MGDoFAccessor<2, 2>;

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
template class MGDoFAccessor<1, 3>;
template class MGDoFAccessor<2, 3>;
template class MGDoFAccessor<3, 3>;

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
