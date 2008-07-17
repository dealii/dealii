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
#ifndef __deal2__tria_accessor_templates_h
#define __deal2__tria_accessor_templates_h


#include <base/config.h>
#include <grid/tria.h>
#include <grid/tria_levels.h>
#include <grid/tria_faces.h>
#include <grid/tria_iterator.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.templates.h>
#include <base/geometry_info.h>

#include <cmath>

DEAL_II_NAMESPACE_OPEN

/*------------------------ Functions: TriaAccessor ---------------------------*/

template <int structdim, int dim>
inline
TriaAccessor<structdim,dim>::TriaAccessor (const Triangulation<dim> *parent,
					   const int                 level,
					   const int                 index,
					   const AccessorData       *)
                :
                present_index (index),
                tria (parent)
{

				   // non-cells have no level, so a 0
				   // should have been passed, or a -1
				   // for an end-iterator, or -2 for
				   // an invalid (default constructed)
				   // iterator
  Assert ((level == 0) || (level == -1) || (level == -2),
	  ExcInternalError());
}



template <int structdim, int dim>
inline
void
TriaAccessor<structdim,dim>::copy_from (const TriaAccessor<structdim,dim> &a)
{
  present_index = a.present_index;
  tria = a.tria;
}



template <int structdim, int dim>
inline
bool
TriaAccessor<structdim,dim>::operator == (const TriaAccessor<structdim,dim> &a) const
{
  Assert (tria == a.tria, TriaAccessorExceptions::ExcCantCompareIterators());
  return (present_index == a.present_index);
}



template <int structdim, int dim>
inline
bool
TriaAccessor<structdim,dim>::operator != (const TriaAccessor<structdim,dim> &a) const
{
  Assert (tria == a.tria, TriaAccessorExceptions::ExcCantCompareIterators());
  return (present_index != a.present_index);
}



template <int structdim, int dim>
inline
int
TriaAccessor<structdim,dim>::level ()
{
  return 0;
}



template <int structdim, int dim>
inline
int
TriaAccessor<structdim,dim>::index () const
{
  return present_index;
}



template <int structdim, int dim>
inline
IteratorState::IteratorStates
TriaAccessor<structdim,dim>::state () const
{
  if (present_index>=0)
    return IteratorState::valid;
  else
    if (present_index==-1)
      return IteratorState::past_the_end;
    else
      return IteratorState::invalid;
}



template <int structdim, int dim>
inline
const Triangulation<dim> &
TriaAccessor<structdim,dim>::get_triangulation () const
{
  return *tria;
}



/*------------------------ Functions: TriaAccessor<dim,dim> ---------------------------*/

template <int dim>
inline
TriaAccessor<dim,dim>::TriaAccessor (const Triangulation<dim> *parent,
				     const int                 level,
				     const int                 index,
				     const AccessorData       *)
                :
                present_level (level),
                present_index (index),
                tria (parent)
{}



template <int dim>
inline
void
TriaAccessor<dim,dim>::copy_from (const TriaAccessor<dim,dim> &a)
{
  present_level = a.present_level;
  present_index = a.present_index;
  tria = a.tria;
}



template <int dim>
inline
bool
TriaAccessor<dim,dim>::operator == (const TriaAccessor<dim,dim> &a) const
{
  Assert (tria == a.tria, TriaAccessorExceptions::ExcCantCompareIterators());
  return ((present_index == a.present_index) &&
	  (present_level == a.present_level));
}



template <int dim>
inline
bool
TriaAccessor<dim,dim>::operator != (const TriaAccessor<dim,dim> &a) const
{
  Assert (tria == a.tria, TriaAccessorExceptions::ExcCantCompareIterators());
  return ((present_index != a.present_index) ||
	  (present_level != a.present_level));
}



template <int dim>
inline
int
TriaAccessor<dim,dim>::level () const
{
  return present_level;
}



template <int dim>
inline
int
TriaAccessor<dim,dim>::index () const
{
  return present_index;
}



template <int dim>
inline
IteratorState::IteratorStates
TriaAccessor<dim,dim>::state () const
{
  if ((present_level>=0) && (present_index>=0))
    return IteratorState::valid;
  else
    if ((present_index==-1) && (present_index==-1))
      return IteratorState::past_the_end;
    else
      return IteratorState::invalid;
}



template <int dim>
inline
const Triangulation<dim> &
TriaAccessor<dim,dim>::get_triangulation () const
{
  return *tria;
}


/*------------------------ Functions: LineAccessor ---------------------------*/


template <int dim>
inline
TriaObjectAccessor<1,dim>::TriaObjectAccessor (
  const Triangulation<dim> *parent,
  const int                 level,
  const int                 index,
  const AccessorData       *local_data)
                :
                TriaAccessor<1,dim> (parent, level, index, local_data)
{
  if (dim!=1)
    Assert(level <= 0, ExcInternalError());
}


// first implement the lines functions as they are used by many
// subsequent accessor functions
template <>
inline
internal::Triangulation::TriaObjects<internal::Triangulation::TriaObject<1> > &
TriaObjectAccessor<1,1>::objects() const
{
  return this->tria->levels[this->present_level]->cells;
}



template <int dim>
inline
internal::Triangulation::TriaObjects<internal::Triangulation::TriaObject<1> > &
TriaObjectAccessor<1,dim>::objects() const
{
  return this->tria->faces->lines;
}



template <int dim>
inline
bool
TriaObjectAccessor<1,dim>::used () const
{
  Assert (this->state() == IteratorState::valid,
	  TriaAccessorExceptions::ExcDereferenceInvalidObject());
  return objects().used[this->present_index];
}



template <int dim>
inline
bool
TriaObjectAccessor<1,dim>::user_flag_set () const
{
  Assert (this->used(), TriaAccessorExceptions::ExcCellNotUsed());
  return objects().user_flags[this->present_index];
}



template <int dim>
inline
void
TriaObjectAccessor<1,dim>::set_user_flag () const 
{
  Assert (this->used(), TriaAccessorExceptions::ExcCellNotUsed());
  objects().user_flags[this->present_index] = true;
}



template <int dim>
inline
void
TriaObjectAccessor<1,dim>::clear_user_flag () const
{
  Assert (this->used(), TriaAccessorExceptions::ExcCellNotUsed());
  objects().user_flags[this->present_index] = false;
}


template <int dim>
inline
RefinementCase<1>
TriaObjectAccessor<1, dim>::refinement_case() const
{
  Assert (this->state() == IteratorState::valid,
	  TriaAccessorExceptions::ExcDereferenceInvalidObject());

  return objects().children[this->present_index] != -1 ?
    RefinementCase<1>::cut_x : RefinementCase<1>::no_refinement;
}



template <int dim>
inline
bool
TriaObjectAccessor<1,dim>::has_children () const
{
  Assert (this->state() == IteratorState::valid,
	  TriaAccessorExceptions::ExcDereferenceInvalidObject());
  return (objects().children[this->present_index] != -1);
}



template <int dim>
inline
int
TriaObjectAccessor<1,dim>::child_index (const unsigned int i) const
{
  Assert (i<2, ExcIndexRange(i,0,2));
  Assert (has_children(), TriaAccessorExceptions::ExcCellHasNoChildren());
  return objects().children[this->present_index]+i;
}



template <int dim>
inline
int
TriaObjectAccessor<1,dim>::isotropic_child_index (const unsigned int i) const
{
  return child_index(i);
}



template <int dim>
inline
TriaIterator<dim,TriaObjectAccessor<1,dim> >
TriaObjectAccessor<1,dim>::child (const unsigned int i) const
{
  Assert (i<2, ExcIndexRange(i,0,2));

  TriaIterator<dim,TriaObjectAccessor<1,dim> >
    q (this->tria,
       (dim == 1 ? this->level() + 1 : 0),
       child_index (i));

  Assert ((q.state() == IteratorState::past_the_end) || q->used(),
	  TriaAccessorExceptions::ExcUnusedCellAsChild());

  return q;
}



template <int dim>
inline
TriaIterator<dim,TriaObjectAccessor<1,dim> >
TriaObjectAccessor<1,dim>::isotropic_child (const unsigned int i) const
{
				   // no anisotropic refinement in 1D
  return child(i);
}



template <int dim>
inline
unsigned int
TriaObjectAccessor<1,dim>::n_children () const
{
  return (has_children() ? GeometryInfo<1>::max_children_per_cell : 0);
}



template <int dim>
inline
unsigned int
TriaObjectAccessor<1,dim>::middle_vertex_index () const
{
  if (has_children())
    return child(0)->vertex_index(1);
  return numbers::invalid_unsigned_int;
}



template <int dim>
inline
unsigned int
TriaObjectAccessor<1,dim>::max_refinement_depth () const
{
  if (!has_children())
    return 0;

  const unsigned int depths[2] = { child(0)->max_refinement_depth() + 1,
				   child(1)->max_refinement_depth() + 1  };
  return std::max (depths[0], depths[1]);
}



template <int dim>
inline
bool
TriaObjectAccessor<1,dim>::face_orientation (const unsigned int) const
{
  return true;
}


template <int dim>
inline
bool
TriaObjectAccessor<1,dim>::face_flip (const unsigned int) const
{
  return false;
}


template <int dim>
inline
bool
TriaObjectAccessor<1,dim>::face_rotation (const unsigned int) const
{
  return false;
}


template <int dim>
inline
bool
TriaObjectAccessor<1,dim>::line_orientation (const unsigned int) const
{
  return true;
}



template <int dim>
inline
void
TriaObjectAccessor<1,dim>::operator ++ ()
{
  ++this->present_index;
				   // is index still in the range of
				   // the vector? (note that we don't
				   // have to set the level, since
				   // dim!=1 and the object therefore
				   // has no level)
  if (this->present_index
      >=
      static_cast<int>(objects().cells.size()))
    this->present_index = -1;
}



template <>
inline
void
TriaObjectAccessor<1,1>::operator ++ ()
{
  ++this->present_index;
				   // is index still in the range of
				   // the vector?
  while (this->present_index
	 >=
	 static_cast<int>(this->tria->levels[this->present_level]->cells.cells.size()))
    {
				       // no -> go one level up until we find
				       // one with more than zero cells
      ++this->present_level;
      this->present_index = 0;
				       // highest level reached?
      if (this->present_level >= static_cast<int>(this->tria->levels.size()))
	{
					   // return with past the end pointer
	  this->present_level = this->present_index = -1;
	  return;
	}
    }
}



template <int dim>
inline
void
TriaObjectAccessor<1,dim>::operator -- ()
{
  --this->present_index;
				   // is index still in the range of
				   // the vector? (note that we don't
				   // have to set the level, since
				   // dim!=1 and the object therefore
				   // has no level)
  if (this->present_index < 0) 
    this->present_index = -1;
  return;
}


template <>
inline
void
TriaObjectAccessor<1,1>::operator -- ()
{
  --this->present_index;
				   // is index still in the range of
				   // the vector?
  while (this->present_index < 0) 
    {
				       // no -> go one level down
      --this->present_level;
				       // lowest level reached?
      if (this->present_level == -1) 
	{
					   // return with past the end pointer
	  this->present_level = this->present_index = -1;
	  return;
	}
					 // else
      this->present_index = this->tria->levels[this->present_level]->cells.cells.size()-1;
    }
}


/*------------------------ Functions: QuadAccessor ---------------------------*/


template <int dim>
inline
TriaObjectAccessor<2,dim>::TriaObjectAccessor (
  const Triangulation<dim> *parent,
  const int                 level,
  const int                 index,
  const AccessorData       *local_data)
                :
                TriaAccessor<2,dim> (parent, level, index, local_data)
{
  if (dim!=2)
    Assert(level <= 0, ExcInternalError());
}



// first implement the quads functions as they are used by many
// subsequent accessor functions
template <>
inline
internal::Triangulation::TriaObjects<internal::Triangulation::TriaObject<2> > &
TriaObjectAccessor<2,2>::objects() const
{
  return this->tria->levels[this->present_level]->cells;
}



template <int dim>
inline
internal::Triangulation::TriaObjects<internal::Triangulation::TriaObject<2> > &
TriaObjectAccessor<2,dim>::objects() const
{
  return this->tria->faces->quads;
}



template <int dim>
inline
bool
TriaObjectAccessor<2,dim>::used () const
{
  Assert (this->state() == IteratorState::valid,
	  TriaAccessorExceptions::ExcDereferenceInvalidObject());
  return objects().used[this->present_index];
}



template <int dim>
inline
bool
TriaObjectAccessor<2,dim>::user_flag_set () const
{
  Assert (this->used(), TriaAccessorExceptions::ExcCellNotUsed());
  return objects().user_flags[this->present_index];
}



template <int dim>
inline
void
TriaObjectAccessor<2,dim>::set_user_flag () const
{
  Assert (this->used(), TriaAccessorExceptions::ExcCellNotUsed());
  objects().user_flags[this->present_index] = true;
}



template <int dim>
inline
void
TriaObjectAccessor<2,dim>::clear_user_flag () const
{
  Assert (this->used(), TriaAccessorExceptions::ExcCellNotUsed());
  objects().user_flags[this->present_index] = false;
}



template <int dim>
inline
TriaIterator<dim,TriaObjectAccessor<1,dim> >
TriaObjectAccessor<2,dim>::line (const unsigned int i) const
{
  return
    TriaIterator<dim,TriaObjectAccessor<1,dim> >
    (
      this->tria,
      0,
      line_index (i)
    );
}



template <int dim>
inline
unsigned int
TriaObjectAccessor<2,dim>::line_index (const unsigned int i) const
{
  Assert (i<4, ExcIndexRange(i,0,4));
  return objects().cells[this->present_index].face(i);
}



template <int dim>
inline
RefinementCase<2>
TriaObjectAccessor<2, dim>::refinement_case () const
{
  Assert (this->state() == IteratorState::valid,
	  TriaAccessorExceptions::ExcDereferenceInvalidObject());
  Assert (static_cast<unsigned int> (this->present_index) <
	  objects().refinement_cases.size(),
	  ExcIndexRange(this->present_index, 0,
			objects().refinement_cases.size()));

  return objects().refinement_cases[this->present_index];
}



template <int dim>
inline
void
TriaObjectAccessor<2, dim>::set_refinement_case (const RefinementCase<2> &refinement_case) const
{
  Assert (this->state() == IteratorState::valid,
	  TriaAccessorExceptions::ExcDereferenceInvalidObject());
  Assert (static_cast<unsigned int> (this->present_index) <
	  objects().refinement_cases.size(),
	  ExcIndexRange(this->present_index, 0,
			objects().refinement_cases.size()));

  objects().refinement_cases[this->present_index] = refinement_case;
}


template <int dim>
inline
void
TriaObjectAccessor<2, dim>::clear_refinement_case () const
{
  Assert (this->state() == IteratorState::valid,
	  TriaAccessorExceptions::ExcDereferenceInvalidObject());
  Assert (static_cast<unsigned int> (this->present_index) <
	  objects().refinement_cases.size(),
	  ExcIndexRange(this->present_index, 0,
			objects().refinement_cases.size()));

  objects().refinement_cases[this->present_index] = RefinementCase<2>::no_refinement;
}


template <int dim>
inline
bool
TriaObjectAccessor<2,dim>::has_children () const
{
  Assert (this->state() == IteratorState::valid,
	  TriaAccessorExceptions::ExcDereferenceInvalidObject());
  return (objects().children[2*this->present_index] != -1);
}



template <int dim>
inline
int TriaObjectAccessor<2,dim>::child_index (const unsigned int i) const
{
  Assert (i<4, ExcIndexRange(i,0,4));
  Assert (has_children(), TriaAccessorExceptions::ExcCellHasNoChildren());
  return objects().children[2*this->present_index+i/2]+i%2;
}



template <int dim>
inline
int TriaObjectAccessor<2,dim>::isotropic_child_index (const unsigned int i) const
{
  switch (static_cast<unsigned char> (refinement_case()))
    {
      case RefinementCase<dim>::cut_x:
					     // this cell is refined with cut_x,
					     // so the child has to be refined
					     // with cut_y
 	    if(child(i%2)->refinement_case()==RefinementCase<dim>::cut_y)
	      return child(i%2)->child_index(i/2);
	    else
	     Assert(false, ExcMessage("This cell has no grandchildren equivalent to isotropic refinement"));
	    break;
      case RefinementCase<dim>::cut_y:
					     // this cell is refined with cut_y,
					     // so the child has to be refined
					     // with cut_x
 	    if (child(i/2)->refinement_case()==RefinementCase<dim>::cut_x)
	      return child(i/2)->child_index(i%2);
	    else
	      Assert(false, ExcMessage("This cell has no grandchildren equivalent to isotropic refinement"));
	    break;
      case RefinementCase<dim>::cut_xy:
	    return child_index(i);
	    break;
      default:
	    Assert(false, TriaAccessorExceptions::ExcCellHasNoChildren());
	    break;
    }
  return -1;
}



template <int dim>
inline
TriaIterator<dim,TriaObjectAccessor<2,dim> >
TriaObjectAccessor<2,dim>::child (const unsigned int i) const
{
  Assert (i<4, ExcIndexRange(i,0,4));
  
  TriaIterator<dim,TriaObjectAccessor<2,dim> >
    q (this->tria,
       (dim == 2 ? this->level() + 1 : 0),
       child_index (i));
  
  Assert ((q.state() == IteratorState::past_the_end) || q->used(),
	  TriaAccessorExceptions::ExcUnusedCellAsChild());

  return q;
}



template <int dim>
inline
TriaIterator<dim,TriaObjectAccessor<2,dim> >
TriaObjectAccessor<2,dim>::isotropic_child (const unsigned int i) const
{
  switch (static_cast<unsigned char> (refinement_case()))
    {
      case RefinementCase<dim>::cut_x:
					     // this cell is refined with cut_x,
					     // so the child has to be refined
					     // with cut_y
 	    Assert(child(i%2)->refinement_case()==RefinementCase<dim>::cut_y,
 		   ExcMessage("This cell has no grandchildren equivalent to isotropic refinement"));
	    return child(i%2)->child(i/2);
	    break;
      case RefinementCase<dim>::cut_y:
					     // this cell is refined with cut_y,
					     // so the child has to be refined
					     // with cut_x
 	    Assert(child(i/2)->refinement_case()==RefinementCase<dim>::cut_x,
 		   ExcMessage("This cell has no grandchildren equivalent to isotropic refinement"));
	    return child(i/2)->child(i%2);
	    break;
      default:
 	    Assert(refinement_case()==RefinementCase<dim>::cut_xy,
 		   TriaAccessorExceptions::ExcCellHasNoChildren());
	    break;
    }
  return child(i);
}



template <int dim>
inline
unsigned int
TriaObjectAccessor<2,dim>::n_children () const
{
  Assert (static_cast<unsigned int> (this->present_index) <
	  objects().refinement_cases.size(),
	  ExcIndexRange(this->present_index, 0,
			objects().refinement_cases.size()));
  return GeometryInfo<2>::n_children(refinement_case());
}



template <int dim>
inline
unsigned int
TriaObjectAccessor<2,dim>::middle_vertex_index () const
{
  switch (static_cast<unsigned char> (refinement_case()))
    {
      case RefinementCase<dim>::cut_x:
	    return child(0)->line(1)->middle_vertex_index();
	    break;
      case RefinementCase<dim>::cut_y:
	    return child(0)->line(3)->middle_vertex_index();
	    break;
      case RefinementCase<dim>::cut_xy:
	    return child(0)->vertex_index(3);
	    break;
      default:
	    break;
    }
  return numbers::invalid_unsigned_int;
}



template <int dim>
inline
unsigned int
TriaObjectAccessor<2,dim>::max_refinement_depth () const
{
  if (!has_children())
    return 0;

  const unsigned int depths[4] = { child(0)->max_refinement_depth() + 1,
				   child(1)->max_refinement_depth() + 1,
				   child(2)->max_refinement_depth() + 1,
				   child(3)->max_refinement_depth() + 1 };
  return std::max (std::max (depths[0], depths[1]),
		   std::max (depths[2], depths[3]));
}



template <int dim>
inline
bool
TriaObjectAccessor<2,dim>::face_orientation (const unsigned int) const
{
  return true;
}


template <int dim>
inline
bool
TriaObjectAccessor<2,dim>::face_flip (const unsigned int) const
{
  return false;
}


template <int dim>
inline
bool
TriaObjectAccessor<2,dim>::face_rotation (const unsigned int) const
{
  return false;
}


template <>
inline
bool
TriaObjectAccessor<2,3>::line_orientation (const unsigned int line) const
{
 				   // we cannot use the objects() function here,
				   // since it returns a reference to
				   // TriaObjects<Quad>, but we need a
				   // (reference to) TriaObjectsQuad3D
  Assert (used(), TriaAccessorExceptions::ExcCellNotUsed());
  return this->tria->faces->quads.face_orientation(this->present_index, line);
}


template <int dim>
inline
bool
TriaObjectAccessor<2,dim>::line_orientation (const unsigned int) const
{
  return true;
}


template <int dim>
inline
void
TriaObjectAccessor<2,dim>::operator ++ ()
{
  ++this->present_index;
				   // is index still in the range of
				   // the vector? (note that we don't
				   // have to set the level, since
				   // dim!=1 and the object therefore
				   // has no level)
  if (this->present_index
      >=
      static_cast<int>(objects().cells.size()))
    this->present_index = -1;
}



template <>
inline
void
TriaObjectAccessor<2,2>::operator ++ ()
{
  ++this->present_index;
				   // is index still in the range of
				   // the vector?
  while (this->present_index
	 >=
	 static_cast<int>(objects().cells.size()))
    {
				       // no -> go one level up
      ++this->present_level;
      this->present_index = 0;
				       // highest level reached?
      if (this->present_level >= static_cast<int>(this->tria->levels.size()))
	{
					   // return with past the end pointer
	  this->present_level = this->present_index = -1;
	  return;
	}
    }
}



template <int dim>
inline
void
TriaObjectAccessor<2,dim>::operator -- ()
{
  --this->present_index;
				   // is index still in the range of
				   // the vector? (note that we don't
				   // have to set the level, since
				   // dim!=1 and the object therefore
				   // has no level)
  if (this->present_index < 0)
    this->present_index = -1;
}


template <>
inline
void
TriaObjectAccessor<2,2>::operator -- ()
{
  --this->present_index;
				   // is index still in the range of
				   // the vector?
  while (this->present_index < 0) 
    {
				       // no -> go one level down
      --this->present_level;
				       // lowest level reached?
      if (this->present_level == -1) 
	{
					   // return with past the end pointer
	  this->present_level = this->present_index = -1;
	  return;
	}
				       // else
      this->present_index = objects().cells.size()-1;
    }
}


/*------------------------ Functions: HexAccessor ---------------------------*/


template <int dim>
inline
TriaObjectAccessor<3,dim>::TriaObjectAccessor (
  const Triangulation<dim> *parent,
  const int                 level,
  const int                 index,
  const AccessorData       *local_data)
                :
                TriaAccessor<3,dim> (parent, level, index, local_data)
{
  if (dim!=3)
    Assert(level <= 0, ExcInternalError());
}



template <>
inline
bool
TriaObjectAccessor<3,3>::used () const
{
  Assert (this->state() == IteratorState::valid,
	  TriaAccessorExceptions::ExcDereferenceInvalidObject());
   return this->tria->levels[this->present_level]->cells.used[this->present_index];
}



template <>
inline
bool
TriaObjectAccessor<3,3>::user_flag_set () const
{
  Assert (this->used(), TriaAccessorExceptions::ExcCellNotUsed());
  return this->tria->levels[this->present_level]->cells.user_flags[this->present_index];
}



template <>
inline
void
TriaObjectAccessor<3,3>::set_user_flag () const
{
  Assert (this->used(), TriaAccessorExceptions::ExcCellNotUsed());
  this->tria->levels[this->present_level]->cells.user_flags[this->present_index] = true;
}



template <>
inline
void TriaObjectAccessor<3,3>::clear_user_flag () const
{
  Assert (this->used(), TriaAccessorExceptions::ExcCellNotUsed());
  this->tria->levels[this->present_level]->cells.user_flags[this->present_index] = false;
}



template <int dim>
inline
TriaIterator<dim,TriaObjectAccessor<1,dim> >
TriaObjectAccessor<3,dim>::line (const unsigned int i) const
{
  Assert (this->used(), TriaAccessorExceptions::ExcCellNotUsed());
  return
    TriaIterator<dim,TriaObjectAccessor<1,dim> >
    (
      this->tria,
      0,
      line_index (i)
    );
}



template <int dim>
inline
TriaIterator<dim,TriaObjectAccessor<2,dim> >
TriaObjectAccessor<3,dim>::quad (const unsigned int i) const
{
  Assert (this->used(), TriaAccessorExceptions::ExcCellNotUsed());
  return
    TriaIterator<dim,TriaObjectAccessor<2,dim> >
    (
      this->tria,
      0,
      quad_index (i)
    );
}



template <int dim>
inline
unsigned int
TriaObjectAccessor<3,dim>::line_index (const unsigned int i) const
{
  Assert (i<12, ExcIndexRange(i,0,12));

                                   // get the line index by asking the
                                   // quads. first assume standard orientation
                                   //
                                   // set up a table that for each
                                   // line describes a) from which
                                   // quad to take it, b) which line
                                   // therein it is if the face is
                                   // oriented correctly
  static const unsigned int lookup_table[12][2] =
    { { 4, 0 }, // take first four lines from bottom face
      { 4, 1 },
      { 4, 2 },
      { 4, 3 },

      { 5, 0 }, // second four lines from top face
      { 5, 1 },
      { 5, 2 },
      { 5, 3 },

      { 0, 0 }, // the rest randomly
      { 1, 0 },
      { 0, 1 },
      { 1, 1 }};

				   // respect non-standard faces by calling the
				   // reordering function from GeometryInfo

  const unsigned int quad_index=lookup_table[i][0];
  const unsigned int std_line_index=lookup_table[i][1];
  
  const unsigned int line_index=GeometryInfo<dim>::standard_to_real_face_line(
    std_line_index,
    face_orientation(quad_index),
    face_flip(quad_index),
    face_rotation(quad_index));

  return (this->quad(quad_index)
          ->line_index(line_index));
}



template <>
inline
unsigned int
TriaObjectAccessor<3,3>::quad_index (const unsigned int i) const
{
  Assert (i<6, ExcIndexRange(i,0,6));

  return this->tria->levels[this->present_level]->cells.cells[this->present_index].face(i);
}



template <>
inline
RefinementCase<3>
TriaObjectAccessor<3, 3>::refinement_case () const
{
  Assert (this->state() == IteratorState::valid,
	  TriaAccessorExceptions::ExcDereferenceInvalidObject());
  Assert (static_cast<unsigned int> (this->present_index) <
	  this->tria->levels[this->present_level]->cells.refinement_cases.size(),
	  ExcIndexRange(this->present_index, 0,
			this->tria->levels[this->present_level]->
			cells.refinement_cases.size()));

  return this->tria->levels[this->present_level]->cells.refinement_cases[this->present_index];
}



template <>
inline
void
TriaObjectAccessor<3, 3>::set_refinement_case (const RefinementCase<3> &refinement_case) const
{
  Assert (static_cast<unsigned int> (this->present_index) <
	  this->tria->levels[this->present_level]->cells.refinement_cases.size(),
	  ExcIndexRange(this->present_index, 0,
			this->tria->levels[this->present_level]->
			cells.refinement_cases.size()));

  this->tria->levels[this->present_level]->
    cells.refinement_cases[this->present_index] = refinement_case;
}


template <>
inline
void
TriaObjectAccessor<3, 3>::clear_refinement_case () const
{
  Assert (this->state() == IteratorState::valid,
	  TriaAccessorExceptions::ExcDereferenceInvalidObject());
  Assert (static_cast<unsigned int> (this->present_index) <
	  this->tria->levels[this->present_level]->cells.refinement_cases.size(),
	  ExcIndexRange(this->present_index, 0,
			this->tria->levels[this->present_level]->
			cells.refinement_cases.size()));

  this->tria->levels[this->present_level]->
    cells.refinement_cases[this->present_index] = RefinementCase<3>::no_refinement;
}



template<>
inline
bool
TriaObjectAccessor<3,3>::has_children () const
{
  Assert (this->state() == IteratorState::valid,
	  TriaAccessorExceptions::ExcDereferenceInvalidObject());
    return (this->tria->levels[this->present_level]->cells.children[4*this->present_index] != -1);
}


template <>
inline
int TriaObjectAccessor<3,3>::child_index (const unsigned int i) const
{
  Assert (i<8, ExcIndexRange(i,0,8));
  Assert (has_children(), TriaAccessorExceptions::ExcCellHasNoChildren());
  return this->tria->levels[this->present_level]->cells.children[4*this->present_index+i/2]+i%2;
}



template <>
inline
int TriaObjectAccessor<3,3>::isotropic_child_index (const unsigned int i) const
{
  AssertThrow(false, ExcNotImplemented());
  return child_index(i);
}



template <>
inline
TriaIterator<3,TriaObjectAccessor<3,3> >
TriaObjectAccessor<3,3>::child (const unsigned int i) const
{
  const int dim=3;
  Assert (i<8, ExcIndexRange(i,0,8));
  
  TriaIterator<dim,TriaObjectAccessor<3,dim> >
    q (this->tria, this->present_level+1, child_index (i));
  
  Assert ((q.state() == IteratorState::past_the_end) || q->used(),
	  TriaAccessorExceptions::ExcUnusedCellAsChild());

  return q;
}



template <>
inline
TriaIterator<3,TriaObjectAccessor<3,3> >
TriaObjectAccessor<3,3>::isotropic_child (const unsigned int i) const
{
  AssertThrow(false, ExcNotImplemented());
  return child(i);
}



template <>
inline
unsigned int
TriaObjectAccessor<3,3>::n_children () const
{
  Assert (static_cast<unsigned int> (this->present_index) <
	  this->tria->levels[this->present_level]->cells.refinement_cases.size(),
	  ExcIndexRange(this->present_index, 0,
			this->tria->levels[this->present_level]->
			cells.refinement_cases.size()));

  return GeometryInfo<3>::n_children(refinement_case());
}



template <int dim>
inline
unsigned int
TriaObjectAccessor<3,dim>::middle_vertex_index () const
{
  switch (static_cast<unsigned char> (refinement_case()))
    {
      case RefinementCase<dim>::cut_x:
	    return child(0)->quad(1)->middle_vertex_index();
	    break;
      case RefinementCase<dim>::cut_y:
	    return child(0)->quad(3)->middle_vertex_index();
	    break;
      case RefinementCase<dim>::cut_z:
	    return child(0)->quad(5)->middle_vertex_index();
	    break;
      case RefinementCase<dim>::cut_xy:
	    return child(0)->line(11)->middle_vertex_index();
	    break;
      case RefinementCase<dim>::cut_xz:
	    return child(0)->line(5)->middle_vertex_index();
	    break;
      case RefinementCase<dim>::cut_yz:
	    return child(0)->line(7)->middle_vertex_index();
	    break;
      case RefinementCase<dim>::cut_xyz:
	    return child(0)->vertex_index(7);
	    break;
      default:
	    break;
    }
  return numbers::invalid_unsigned_int;
}



template <int dim>
inline
unsigned int
TriaObjectAccessor<3,dim>::max_refinement_depth () const
{
  if (!has_children())
    return 0;

  const unsigned int depths[8] = { child(0)->max_refinement_depth() + 1,
				   child(1)->max_refinement_depth() + 1,
				   child(2)->max_refinement_depth() + 1,
				   child(3)->max_refinement_depth() + 1,
				   child(4)->max_refinement_depth() + 1,
				   child(5)->max_refinement_depth() + 1,
				   child(6)->max_refinement_depth() + 1,
				   child(7)->max_refinement_depth() + 1  };
  return std::max (std::max (std::max (depths[0], depths[1]),
			     std::max (depths[2], depths[3])),
		   std::max (std::max (depths[4], depths[5]),
			     std::max (depths[6], depths[7])));
}



template <>
inline
bool
TriaObjectAccessor<3, 3>::face_orientation (const unsigned int face) const
{
  Assert (used(), TriaAccessorExceptions::ExcCellNotUsed());  
  return (this->tria->levels[this->present_level]
	  ->cells.face_orientation(this->present_index, face));
}


template <>
inline
bool
TriaObjectAccessor<3, 3>::face_flip (const unsigned int face) const
{
  Assert (used(), TriaAccessorExceptions::ExcCellNotUsed());
  Assert (face<GeometryInfo<3>::faces_per_cell,
          ExcIndexRange (face, 0, GeometryInfo<3>::faces_per_cell));
  Assert (this->present_index * GeometryInfo<3>::faces_per_cell + face
	  < this->tria->levels[this->present_level]
	  ->cells.face_flips.size(),
	  ExcInternalError());
  
  return (this->tria->levels[this->present_level]
	  ->cells.face_flips[this->present_index *
			     GeometryInfo<3>::faces_per_cell
			     + face]);
}


template <>
inline
bool
TriaObjectAccessor<3, 3>::face_rotation (const unsigned int face) const
{
  Assert (used(), TriaAccessorExceptions::ExcCellNotUsed());
  Assert (face<GeometryInfo<3>::faces_per_cell,
          ExcIndexRange (face, 0, GeometryInfo<3>::faces_per_cell));
  Assert (this->present_index * GeometryInfo<3>::faces_per_cell + face
	  < this->tria->levels[this->present_level]
	  ->cells.face_rotations.size(),
	  ExcInternalError());
      
  return (this->tria->levels[this->present_level]
	  ->cells.face_rotations[this->present_index *
				 GeometryInfo<3>::faces_per_cell
				 + face]);
}


template <>
inline
bool
TriaObjectAccessor<3, 3>::line_orientation (const unsigned int line) const
{
  const int dim=3;
  Assert (used(), TriaAccessorExceptions::ExcCellNotUsed());
  Assert (line<GeometryInfo<3>::lines_per_cell,
          ExcIndexRange (line, 0, GeometryInfo<3>::lines_per_cell));

                                   // get the line index by asking the
                                   // quads. first assume standard orientation
                                   //
                                   // set up a table that for each
                                   // line describes a) from which
                                   // quad to take it, b) which line
                                   // therein it is if the face is
                                   // oriented correctly
  static const unsigned int lookup_table[12][2] =
    { { 4, 0 }, // take first four lines from bottom face
      { 4, 1 },
      { 4, 2 },
      { 4, 3 },

      { 5, 0 }, // second four lines from top face
      { 5, 1 },
      { 5, 2 },
      { 5, 3 },

      { 0, 0 }, // the rest randomly
      { 1, 0 },
      { 0, 1 },
      { 1, 1 }};

  const unsigned int quad_index=lookup_table[line][0];
  const unsigned int std_line_index=lookup_table[line][1];
  
  const unsigned int line_index=GeometryInfo<dim>::standard_to_real_face_line(
    std_line_index,
    face_orientation(quad_index),
    face_flip(quad_index),
    face_rotation(quad_index));
  
				   // now we got to the correct line and ask
				   // the quad for its line_orientation. however, if
				   // the face is rotated, it might be possible,
				   // that a standard orientation of the line
				   // with respect to the face corrsponds to a
				   // non-standard orientation for the line with
				   // respect to the cell.
				   //
				   // set up a table indicating if the two
				   // standard orientations coincide
				   //
				   // first index: two pairs of lines 0(lines
				   // 0/1) and 1(lines 2/3)
				   //
				   // second index: face_orientation; 0:
				   // opposite normal, 1: standard
				   //
				   // third index: face_flip; 0: standard, 1:
				   // face rotated by 180 degrees
				   //
				   // forth index: face_rotation: 0: standard,
				   // 1: face rotated by 90 degrees

  static const bool bool_table[2][2][2][2] =
    { { { { true, false },   // lines 0/1, face_orientation=false, face_flip=false, face_rotation=false and true
	  { false, true }},  // lines 0/1, face_orientation=false, face_flip=true, face_rotation=false and true
	{ { true, true },    // lines 0/1, face_orientation=true, face_flip=false, face_rotation=false and true
	  { false, false }}},// linea 0/1, face_orientation=true, face_flip=true, face_rotation=false and true

      { { { true, true },    // lines 2/3 ...
	  { false, false }},
	{ { true, false },
	  { false, true }}}};
  
  
  return (this->quad(quad_index)
          ->line_orientation(line_index)
	  == bool_table[std_line_index/2]
	  [face_orientation(quad_index)]
	  [face_flip(quad_index)]
	  [face_rotation(quad_index)]);
}



template <>
inline
void
TriaObjectAccessor<3,3>::operator ++ ()
{
  ++this->present_index;
				   // is index still in the range of
				   // the vector?
    while (this->present_index
	   >=
	   static_cast<int>(this->tria->levels[this->present_level]->cells.cells.size()))
      {
					 // no -> go one level up
	++this->present_level;
	this->present_index = 0;
					 // highest level reached?
	if (this->present_level >= static_cast<int>(this->tria->levels.size()))
	  {
					     // return with past the end pointer
	    this->present_level = this->present_index = -1;
	    return;
	  }
      }
}



template <>
inline
void
TriaObjectAccessor<3,3>::operator -- ()
{
  --this->present_index;
				   // is index still in the range of
				   // the vector?
    while (this->present_index < 0) 
      {
					 // no -> go one level down
	--this->present_level;
					 // lowest level reached?
	if (this->present_level == -1) 
	  {
					     // return with past the end pointer
	    this->present_level = this->present_index = -1;
	    return;
	  }
					 // else
	this->present_index = this->tria->levels[this->present_level]->cells.cells.size()-1;
      }
}


/*------------------------ Functions: CellAccessor<dim> -----------------------*/


template <int dim>
inline
CellAccessor<dim>::CellAccessor (
  const Triangulation<dim> *parent,
  const int                 level,
  const int                 index,
  const AccessorData       *local_data)
                :
                TriaObjectAccessor<dim,dim> (parent, level, index, local_data)
{}



template <>
inline
TriaIterator<1,TriaObjectAccessor<0, 1> >
CellAccessor<1>::face (const unsigned int) const 
{
  Assert (false, ExcImpossibleInDim(1));
  return TriaIterator<1,TriaObjectAccessor<0, 1> >();
}



template <>
inline
Triangulation<2>::face_iterator
CellAccessor<2>::face (const unsigned int i) const 
{
  return this->line(i);
}



template <>
inline
Triangulation<3>::face_iterator
CellAccessor<3>::face (const unsigned int i) const 
{
  return this->quad(i);
}



template <>
inline
unsigned int
CellAccessor<1>::face_index (const unsigned int) const 
{
  Assert (false, ExcImpossibleInDim(1));
  return numbers::invalid_unsigned_int;
}



template <>
inline
unsigned int
CellAccessor<2>::face_index (const unsigned int i) const 
{
  return this->line_index(i);
}



template <>
inline
unsigned int
CellAccessor<3>::face_index (const unsigned int i) const 
{
  return this->quad_index(i);
}



template <int dim>
inline
int
CellAccessor<dim>::neighbor_index (const unsigned int i) const 
{
  Assert (i<GeometryInfo<dim>::faces_per_cell,
	  TriaAccessorExceptions::ExcInvalidNeighbor(i));
  return this->tria->levels[this->present_level]->
    neighbors[this->present_index*GeometryInfo<dim>::faces_per_cell+i].second;
}



template <int dim>
inline
int
CellAccessor<dim>::neighbor_level (const unsigned int i) const
{
  Assert (i<GeometryInfo<dim>::faces_per_cell,
	  TriaAccessorExceptions::ExcInvalidNeighbor(i));
  return this->tria->levels[this->present_level]->
    neighbors[this->present_index*GeometryInfo<dim>::faces_per_cell+i].first;
}



template <int dim>
inline
RefinementCase<dim>
CellAccessor<dim>::refine_flag_set () const
{
  Assert (this->used(), TriaAccessorExceptions::ExcCellNotUsed());
				   // cells flagged for refinement must be active
				   // (the @p set_refine_flag function checks this,
				   // but activity may change when refinement is
				   // executed and for some reason the refine
				   // flag is not cleared).
  Assert (this->active() ||  !this->tria->levels[this->present_level]->refine_flags[this->present_index],
	  ExcRefineCellNotActive());
  return RefinementCase<dim>(this->tria->levels[this->present_level]->refine_flags[this->present_index]);
}



template <int dim>
inline
void
CellAccessor<dim>::set_refine_flag (const RefinementCase<dim> refinement_case) const
{
  Assert (this->used() && this->active(), ExcRefineCellNotActive());
  Assert (!coarsen_flag_set(),
	  ExcCellFlaggedForCoarsening());
  
  this->tria->levels[this->present_level]->refine_flags[this->present_index] = refinement_case;
}



template <int dim>
inline
void
CellAccessor<dim>::clear_refine_flag () const
{
  Assert (this->used() && this->active(), ExcRefineCellNotActive());
  this->tria->levels[this->present_level]->refine_flags[this->present_index] =
    RefinementCase<dim>::no_refinement;
}



template <int dim>
inline
bool
CellAccessor<dim>::flag_for_face_refinement (const unsigned int face_no,
					     const RefinementCase<dim-1> &face_refinement_case) const
{
  Assert (dim>1, ExcImpossibleInDim(dim));
  Assert (face_no<GeometryInfo<dim>::faces_per_cell,
	  ExcIndexRange(face_no,0,GeometryInfo<dim>::faces_per_cell));
  Assert (face_refinement_case < RefinementCase<dim>::isotropic_refinement+1,
	  ExcIndexRange(face_refinement_case,0,RefinementCase<dim>::isotropic_refinement+1));
			    
                                   // the new refinement case is a combination
				   // of the minimum required one for the given
				   // face refinement and the already existing
				   // flagged refinement case
  RefinementCase<dim> old_ref_case = refine_flag_set();
  RefinementCase<dim>		    
    new_ref_case = (old_ref_case
		    | GeometryInfo<dim>::min_cell_refinement_case_for_face_refinement(face_refinement_case,
										      face_no,
										      this->face_orientation(face_no),
										      this->face_flip(face_no),
										      this->face_rotation(face_no)));
  set_refine_flag(new_ref_case);
				   // return, whether we had to change the
				   // refinement flag
  return new_ref_case != old_ref_case;
}



template <int dim>
inline
bool
CellAccessor<dim>::flag_for_line_refinement (const unsigned int line_no) const
{
  Assert (dim>1, ExcImpossibleInDim(dim));
  Assert (line_no<GeometryInfo<dim>::lines_per_cell,
	  ExcIndexRange(line_no,0,GeometryInfo<dim>::lines_per_cell));

				   // the new refinement case is a combination
				   // of the minimum required one for the given
				   // line refinement and the already existing
				   // flagged refinement case
  RefinementCase<dim> old_ref_case=refine_flag_set(),
		   new_ref_case=old_ref_case
				| GeometryInfo<dim>::min_cell_refinement_case_for_line_refinement(line_no);
  set_refine_flag(new_ref_case);
				   // return, whether we had to change the
				   // refinement flag
  return new_ref_case != old_ref_case;
}



template <>
inline
internal::SubfaceCase<1>
CellAccessor<1>::subface_case(const unsigned int) const
{
  Assert(false, ExcImpossibleInDim(1));
  return internal::SubfaceCase<1>::case_none;
}



template <>
inline
internal::SubfaceCase<2>
CellAccessor<2>::subface_case(const unsigned int face_no) const
{
  Assert(active(), TriaAccessorExceptions::ExcCellNotActive());
  Assert(face_no<GeometryInfo<2>::faces_per_cell,
	 ExcIndexRange(face_no,0,GeometryInfo<2>::faces_per_cell));
  return (face(face_no)->has_children()) ? internal::SubfaceCase<2>::case_x : internal::SubfaceCase<2>::case_none;
}



template <>
inline
internal::SubfaceCase<3>
CellAccessor<3>::subface_case(const unsigned int face_no) const
{
  Assert(active(), TriaAccessorExceptions::ExcCellNotActive());
  Assert(face_no<GeometryInfo<3>::faces_per_cell,
	 ExcIndexRange(face_no,0,GeometryInfo<3>::faces_per_cell));
  switch (static_cast<unsigned char> (face(face_no)->refinement_case()))
    {
      case RefinementCase<3>::no_refinement:
	    return internal::SubfaceCase<3>::case_none;
	    break;
      case RefinementCase<3>::cut_x:
	    if (face(face_no)->child(0)->has_children())
	      {
 		Assert(face(face_no)->child(0)->refinement_case()==RefinementCase<3>::cut_y,
 		       ExcInternalError());
		if (face(face_no)->child(1)->has_children())
		  {
 		    Assert(face(face_no)->child(1)->refinement_case()==RefinementCase<3>::cut_y,
 			   ExcInternalError());
		    return internal::SubfaceCase<3>::case_x1y2y;								 
		  }
		else
		  return internal::SubfaceCase<3>::case_x1y;
	      }
	    else
	      {
		if (face(face_no)->child(1)->has_children())
		  {
  		    Assert(face(face_no)->child(1)->refinement_case()==RefinementCase<3>::cut_y,
  			   ExcInternalError());
		    return internal::SubfaceCase<3>::case_x2y;
		  }
		else
		  return internal::SubfaceCase<3>::case_x;
	      }
	    break;
      case RefinementCase<3>::cut_y:
	    if (face(face_no)->child(0)->has_children())
	      {
 		Assert(face(face_no)->child(0)->refinement_case()==RefinementCase<3>::cut_x,
 		       ExcInternalError());
		if (face(face_no)->child(1)->has_children())
		  {
 		    Assert(face(face_no)->child(1)->refinement_case()==RefinementCase<3>::cut_x,
 			   ExcInternalError());
		    return internal::SubfaceCase<3>::case_y1x2x;								 
		  }
		else
		  return internal::SubfaceCase<3>::case_y1x;
	      }
	    else
	      {
		if (face(face_no)->child(1)->has_children())
		  {
  		    Assert(face(face_no)->child(1)->refinement_case()==RefinementCase<3>::cut_x,
  			   ExcInternalError());
		    return internal::SubfaceCase<3>::case_y2x;
		  }
		else
		  return internal::SubfaceCase<3>::case_y;
	      }
	    break;
      case RefinementCase<3>::cut_xy:
	    return internal::SubfaceCase<3>::case_xy;
	    break;
      default:
	    Assert(false, ExcInternalError());
    }
				   // we should never get here
  return internal::SubfaceCase<3>::case_none;
}



template <int dim>
inline
bool
CellAccessor<dim>::coarsen_flag_set () const
{
  Assert (this->used(), TriaAccessorExceptions::ExcCellNotUsed());
				   // cells flagged for coarsening must be active
				   // (the @p set_refine_flag function checks this,
				   // but activity may change when refinement is
				   // executed and for some reason the refine
				   // flag is not cleared).
  Assert (this->active() ||  !this->tria->levels[this->present_level]->coarsen_flags[this->present_index],
	  ExcRefineCellNotActive());
  return this->tria->levels[this->present_level]->coarsen_flags[this->present_index];
}



template <int dim>
inline
void
CellAccessor<dim>::set_coarsen_flag () const
{
  Assert (this->used() && this->active(), ExcRefineCellNotActive());
  Assert (!refine_flag_set(), ExcCellFlaggedForRefinement());
  
  this->tria->levels[this->present_level]->coarsen_flags[this->present_index] = true;
}



template <int dim>
inline
void
CellAccessor<dim>::clear_coarsen_flag () const
{
  Assert (this->used() && this->active(), ExcRefineCellNotActive());
  this->tria->levels[this->present_level]->coarsen_flags[this->present_index] = false;
}



template <int dim>
inline
TriaIterator<dim,CellAccessor<dim> >
CellAccessor<dim>::neighbor (const unsigned int i) const
{
  TriaIterator<dim,CellAccessor<dim> >
    q (this->tria, neighbor_level (i), neighbor_index (i));

  Assert ((q.state() == IteratorState::past_the_end) || q->used(),
	  TriaAccessorExceptions::ExcUnusedCellAsNeighbor());

  return q;
}



template <int dim>
inline
TriaIterator<dim,CellAccessor<dim> >
CellAccessor<dim>::child (const unsigned int i) const
{
  TriaIterator<dim,CellAccessor<dim> >
    q (this->tria, this->present_level+1, this->child_index (i));

  Assert ((q.state() == IteratorState::past_the_end) || q->used(),
	  TriaAccessorExceptions::ExcUnusedCellAsChild());

  return q;
}



template <int dim>
inline
bool
CellAccessor<dim>::active () const
{
  return !this->has_children();
}



template <int dim>
inline
unsigned int
CellAccessor<dim>::neighbor_face_no (const unsigned int neighbor) const
{
  if (dim==1)
    return neighbor_of_neighbor(neighbor);
  else
    {
      const unsigned int n2=neighbor_of_neighbor_internal(neighbor);
      if (n2!=numbers::invalid_unsigned_int)
					 // return this value as the
					 // neighbor is not coarser
	return n2;
      else	
					 // the neighbor is coarser
	return neighbor_of_coarser_neighbor(neighbor).first;
    }
}

DEAL_II_NAMESPACE_CLOSE

#endif
