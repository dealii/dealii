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
internal::Triangulation::TriaObjects<internal::Triangulation::Line> &
TriaObjectAccessor<1,1>::objects() const
{
  return this->tria->levels[this->present_level]->cells;
}



template <int dim>
inline
internal::Triangulation::TriaObjects<internal::Triangulation::Line> &
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
TriaObjectAccessor<1,dim>::child_index (unsigned const int i) const
{
  Assert (i<2, ExcIndexRange(i,0,2));
  Assert (has_children(), TriaAccessorExceptions::ExcCellHasNoChildren());
  return objects().children[this->present_index]+i;
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
unsigned int
TriaObjectAccessor<1,dim>::n_children () const
{
  Assert (has_children()==true, TriaAccessorExceptions::ExcCellHasNoChildren());
  return GeometryInfo<1>::children_per_cell;
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
      static_cast<int>(this->tria->faces->lines.cells.size()))
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
internal::Triangulation::TriaObjects<internal::Triangulation::Quad> &
TriaObjectAccessor<2,2>::objects() const
{
  return this->tria->levels[this->present_level]->cells;
}



template <int dim>
inline
internal::Triangulation::TriaObjects<internal::Triangulation::Quad> &
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
bool
TriaObjectAccessor<2,dim>::has_children () const
{
  Assert (this->state() == IteratorState::valid,
	  TriaAccessorExceptions::ExcDereferenceInvalidObject());
  return (objects().children[this->present_index] != -1);
}



template <int dim>
inline
int TriaObjectAccessor<2,dim>::child_index (const unsigned int i) const
{
  Assert (i<4, ExcIndexRange(i,0,4));
  Assert (has_children(), TriaAccessorExceptions::ExcCellHasNoChildren());
  return objects().children[this->present_index]+i;
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
unsigned int
TriaObjectAccessor<2,dim>::n_children () const
{
  Assert (has_children()==true, TriaAccessorExceptions::ExcCellHasNoChildren());
  return GeometryInfo<2>::children_per_cell;
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
bool
TriaObjectAccessor<3,3>::has_children () const
{
  Assert (this->state() == IteratorState::valid,
	  TriaAccessorExceptions::ExcDereferenceInvalidObject());
    return (this->tria->levels[this->present_level]->cells.children[this->present_index] != -1);
}


template <>
inline
int TriaObjectAccessor<3,3>::child_index (const unsigned int i) const
{
  Assert (i<8, ExcIndexRange(i,0,8));
  Assert (has_children(), TriaAccessorExceptions::ExcCellHasNoChildren());
    return this->tria->levels[this->present_level]->cells.children[this->present_index]+i;
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



template <int dim>
inline
unsigned int
TriaObjectAccessor<3,dim>::n_children () const
{
  Assert (has_children()==true, TriaAccessorExceptions::ExcCellHasNoChildren());
  return GeometryInfo<3>::children_per_cell;
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
  return deal_II_numbers::invalid_unsigned_int;
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
bool
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
  return this->tria->levels[this->present_level]->refine_flags[this->present_index];
}



template <int dim>
inline
void
CellAccessor<dim>::set_refine_flag () const
{
  Assert (this->used() && this->active(), ExcRefineCellNotActive());
  Assert (!coarsen_flag_set(),
	  ExcCellFlaggedForCoarsening());
  
  this->tria->levels[this->present_level]->refine_flags[this->present_index] = true;
}



template <int dim>
inline
void
CellAccessor<dim>::clear_refine_flag () const
{
  Assert (this->used() && this->active(), ExcRefineCellNotActive());
  this->tria->levels[this->present_level]->refine_flags[this->present_index] = false;
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

DEAL_II_NAMESPACE_CLOSE

#endif
