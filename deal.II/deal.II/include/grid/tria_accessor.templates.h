//----------------------------  tria_accessor.templates.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  tria_accessor.templates.h  ---------------------------
#ifndef __deal2__tria_accessor_templates_h
#define __deal2__tria_accessor_templates_h


#include <grid/tria.h>
#include <grid/tria_levels.h>
#include <grid/tria_iterator.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.templates.h>
#include <grid/geometry_info.h>

#include <cmath>


/*------------------------ Functions: TriaAccessor ---------------------------*/


template <int dim>
inline
void
TriaAccessor<dim>::copy_from (const TriaAccessor<dim> &a) {
  present_level = a.present_level;
  present_index = a.present_index;
  tria = a.tria;
};


template <int dim>
inline
bool
TriaAccessor<dim>::operator == (const TriaAccessor<dim> &a) const {
  Assert (tria == a.tria, ExcCantCompareIterators());
  return ((present_index == a.present_index) &&
	  (present_level == a.present_level));
};


template <int dim>
inline
bool
TriaAccessor<dim>::operator != (const TriaAccessor<dim> &a) const {
  Assert (tria == a.tria, ExcCantCompareIterators());
  return ((present_index != a.present_index) ||
	  (present_level != a.present_level));
};


template <int dim>
inline
int
TriaAccessor<dim>::level () const {
  return present_level;
};


template <int dim>
inline
int
TriaAccessor<dim>::index () const {
  return present_index;
};


template <int dim>
inline
IteratorState
TriaAccessor<dim>::state () const {
  if ((present_level>=0) && (present_index>=0))
    return valid;
  else
    if ((present_level==-1) && (present_index==-1))
      return past_the_end;
    else
      return invalid;
};


template <int dim>
inline
const Triangulation<dim> &
TriaAccessor<dim>::get_triangulation () const
{
  return *tria;
};


/*------------------------ Functions: LineAccessor ---------------------------*/


template <int dim>
inline
bool
TriaObjectAccessor<1,dim>::used () const {
  Assert (state() == valid,
	  typename TriaAccessor<dim>::ExcDereferenceInvalidObject());
  return tria->levels[present_level]->lines.used[present_index];
};


template <int dim>
inline
bool
TriaObjectAccessor<1,dim>::user_flag_set () const {
  Assert (used(), typename TriaAccessor<dim>::ExcCellNotUsed());
  return tria->levels[present_level]->lines.user_flags[present_index];
};


template <int dim>
inline
void
TriaObjectAccessor<1,dim>::set_user_flag () const {
  Assert (used(), typename TriaAccessor<dim>::ExcCellNotUsed());
  tria->levels[present_level]->lines.user_flags[present_index] = true;
};


template <int dim>
inline
void
TriaObjectAccessor<1,dim>::clear_user_flag () const {
  Assert (used(), typename TriaAccessor<dim>::ExcCellNotUsed());
  tria->levels[present_level]->lines.user_flags[present_index] = false;
};


template <int dim>
inline
TriaIterator<dim,TriaObjectAccessor<1,dim> >
TriaObjectAccessor<1,dim>::child (const unsigned int i) const {
  Assert (i<2, ExcIndexRange(i,0,2));
  
  TriaIterator<dim,TriaObjectAccessor<1,dim> > q (tria, present_level+1, child_index (i));
  
#ifdef DEBUG
  if (q.state() != past_the_end)
    Assert (q->used(), typename TriaAccessor<dim>::ExcUnusedCellAsChild());
#endif
  return q;
};


template <int dim>
inline
int
TriaObjectAccessor<1,dim>::child_index (unsigned int i) const {
  Assert (i<2, ExcIndexRange(i,0,2));
  return tria->levels[present_level]->lines.children[present_index]+i;
};


template <int dim>
inline
bool
TriaObjectAccessor<1,dim>::has_children () const {
  Assert (state() == valid, typename TriaAccessor<dim>::ExcDereferenceInvalidObject());
  return (tria->levels[present_level]->lines.children[present_index] != -1);
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
};


template <int dim>
inline
void
TriaObjectAccessor<1,dim>::operator ++ () {
  ++present_index;
				   // is index still in the range of
				   // the vector?
  while (present_index
	 >=
	 static_cast<int>(tria->levels[present_level]->lines.lines.size()))
    {
				       // no -> go one level up until we find
				       // one with more than zero cells
      ++present_level;
      present_index = 0;
				       // highest level reached?
      if (present_level >= static_cast<int>(tria->levels.size()))
	{
					   // return with past the end pointer
	  present_level = present_index = -1;
	  return;
	};
    };
};


template <int dim>
inline
void
TriaObjectAccessor<1,dim>::operator -- () {
  --present_index;
				   // is index still in the range of
				   // the vector?
  while (present_index < 0) 
    {
				       // no -> go one level down
      --present_level;
				       // lowest level reached?
      if (present_level == -1) 
	{
					   // return with past the end pointer
	  present_level = present_index = -1;
	  return;
	};
				       // else
      present_index = tria->levels[present_level]->lines.lines.size()-1;
    };
};


/*------------------------ Functions: QuadAccessor ---------------------------*/


template <int dim>
inline
bool
TriaObjectAccessor<2,dim>::used () const {
  Assert (state() == valid,
	  typename TriaAccessor<dim>::ExcDereferenceInvalidObject());
  return tria->levels[present_level]->quads.used[present_index];
};


template <int dim>
inline
bool
TriaObjectAccessor<2,dim>::user_flag_set () const {
  Assert (used(), typename TriaAccessor<dim>::ExcCellNotUsed());
  return tria->levels[present_level]->quads.user_flags[present_index];
};


template <int dim>
inline
void
TriaObjectAccessor<2,dim>::set_user_flag () const {
  Assert (used(), typename TriaAccessor<dim>::ExcCellNotUsed());
  tria->levels[present_level]->quads.user_flags[present_index] = true;
};


template <int dim>
inline
void
TriaObjectAccessor<2,dim>::clear_user_flag () const {
  Assert (used(), typename TriaAccessor<dim>::ExcCellNotUsed());
  tria->levels[present_level]->quads.user_flags[present_index] = false;
};


template <int dim>
inline
TriaIterator<dim,TriaObjectAccessor<1,dim> >
TriaObjectAccessor<2,dim>::line (const unsigned int i) const {
  return
    TriaIterator<dim,TriaObjectAccessor<1,dim> >
    (
      tria,
      present_level,
      line_index (i)
    );
};


template <int dim>
inline
unsigned int
TriaObjectAccessor<2,dim>::line_index (unsigned int i) const {
  Assert (i<4, ExcIndexRange(i,0,4));

  return tria->levels[present_level]->quads.quads[present_index].line(i);
};


template <int dim>
inline
TriaIterator<dim,TriaObjectAccessor<2,dim> >
TriaObjectAccessor<2,dim>::child (const unsigned int i) const {
  Assert (i<4, ExcIndexRange(i,0,4));
  
  TriaIterator<dim,TriaObjectAccessor<2,dim> > q (tria, present_level+1, child_index (i));
  
#ifdef DEBUG
  if (q.state() != past_the_end)
    Assert (q->used(), typename TriaAccessor<dim>::ExcUnusedCellAsChild());
#endif
  return q;
};


template <int dim>
inline
int TriaObjectAccessor<2,dim>::child_index (unsigned int i) const {
  Assert (i<4, ExcIndexRange(i,0,4));
  return tria->levels[present_level]->quads.children[present_index]+i;
};


template <int dim>
inline
bool
TriaObjectAccessor<2,dim>::has_children () const {
  Assert (state() == valid, typename TriaAccessor<dim>::ExcDereferenceInvalidObject());
  return (tria->levels[present_level]->quads.children[present_index] != -1);
};


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
};


template <int dim>
inline
void
TriaObjectAccessor<2,dim>::operator ++ () {
  ++present_index;
				   // is index still in the range of
				   // the vector?
  while (present_index
	 >=
	 static_cast<int>(tria->levels[present_level]->quads.quads.size()))
    {
				       // no -> go one level up
      ++present_level;
      present_index = 0;
				       // highest level reached?
      if (present_level >= static_cast<int>(tria->levels.size()))
	{
					   // return with past the end pointer
	  present_level = present_index = -1;
	  return;
	};
    };
};


template <int dim>
inline
void
TriaObjectAccessor<2,dim>::operator -- () {
  --present_index;
				   // is index still in the range of
				   // the vector?
  while (present_index < 0) 
    {
				       // no -> go one level down
      --present_level;
				       // lowest level reached?
      if (present_level == -1) 
	{
					   // return with past the end pointer
	  present_level = present_index = -1;
	  return;
	};
				       // else
      present_index = tria->levels[present_level]->quads.quads.size()-1;
    };
};


/*------------------------ Functions: HexAccessor ---------------------------*/


template <int dim>
inline
bool
TriaObjectAccessor<3,dim>::used () const {
  Assert (state() == valid,
	  typename TriaAccessor<dim>::ExcDereferenceInvalidObject());
  return tria->levels[present_level]->hexes.used[present_index];
};


template <int dim>
inline
bool
TriaObjectAccessor<3,dim>::user_flag_set () const {
  Assert (used(), typename TriaAccessor<dim>::ExcCellNotUsed());
  return tria->levels[present_level]->hexes.user_flags[present_index];
};


template <int dim>
inline
void
TriaObjectAccessor<3,dim>::set_user_flag () const {
  Assert (used(), typename TriaAccessor<dim>::ExcCellNotUsed());
  tria->levels[present_level]->hexes.user_flags[present_index] = true;
};


template <int dim>
inline
void TriaObjectAccessor<3,dim>::clear_user_flag () const {
  Assert (used(), typename TriaAccessor<dim>::ExcCellNotUsed());
  tria->levels[present_level]->hexes.user_flags[present_index] = false;
};


template <int dim>
inline
TriaIterator<dim,TriaObjectAccessor<1,dim> >
TriaObjectAccessor<3,dim>::line (const unsigned int i) const {
  Assert (used(), typename TriaAccessor<dim>::ExcCellNotUsed());
  Assert (i<12, ExcIndexRange (i,0,12));

				   // egcs 1.1.2 gets into trouble if we
				   // omit the this-> here, if one tries to
				   // inline this function into another
				   // function where the variable name
				   // quad is also used.. It then complains
				   // that the name look-up for
				   // for-loop-variables has changed. using
				   // this-> as here works around the problem
  if (i<4)
    return this->quad(0)->line(i);
  else
    if (i<8)
      return this->quad(1)->line(i-4);
    else
      switch (i) 
	{
	  case 8:
		return this->quad(2)->line(3);
	  case 9:
		return this->quad(2)->line(1);
	  case 10:
		return this->quad(4)->line(1);
	  case 11:
		return this->quad(4)->line(3);
	};
  Assert (false, ExcIndexRange(i,0,12));
  return TriaIterator<dim,TriaObjectAccessor<1,dim> >(tria, -1, -1, 0);
};


template <int dim>
inline
TriaIterator<dim,TriaObjectAccessor<2,dim> >
TriaObjectAccessor<3,dim>::quad (const unsigned int i) const {
  Assert (used(), typename TriaAccessor<dim>::ExcCellNotUsed());
  return
    TriaIterator<dim,TriaObjectAccessor<2,dim> >
    (
      tria,
      present_level,
      quad_index (i)
    );
};


template <int dim>
inline
unsigned int
TriaObjectAccessor<3,dim>::line_index (unsigned int i) const {
  Assert (i<12, ExcIndexRange(i,0,12));

  if (i<4)
    return quad(0)->line_index(i);
  else
    if (i<8)
      return quad(1)->line_index(i-4);
    else
      switch (i) 
	{
	  case 8:
		return quad(2)->line_index(3);
	  case 9:
		return quad(2)->line_index(1);
	  case 10:
		return quad(4)->line_index(1);
	  case 11:
		return quad(4)->line_index(3);
	};
  Assert (false, ExcIndexRange(i,0,12));
  return 0;
};


template <int dim>
inline
unsigned int
TriaObjectAccessor<3,dim>::quad_index (unsigned int i) const
{
  Assert (i<6, ExcIndexRange(i,0,6));

  return tria->levels[present_level]->hexes.hexes[present_index].quad(i);
};


template <int dim>
inline
TriaIterator<dim,TriaObjectAccessor<3,dim> >
TriaObjectAccessor<3,dim>::child (const unsigned int i) const
{
  Assert (i<8, ExcIndexRange(i,0,8));
  
  TriaIterator<dim,TriaObjectAccessor<3,dim> > q (tria, present_level+1, child_index (i));
  
#ifdef DEBUG
  if (q.state() != past_the_end)
    Assert (q->used(), typename TriaAccessor<dim>::ExcUnusedCellAsChild());
#endif
  return q;
};


template <int dim>
inline
int TriaObjectAccessor<3,dim>::child_index (unsigned int i) const {
  Assert (i<8, ExcIndexRange(i,0,8));
  return tria->levels[present_level]->hexes.children[present_index]+i;
};


template <int dim>
bool TriaObjectAccessor<3,dim>::has_children () const {
  Assert (state() == valid, typename TriaAccessor<dim>::ExcDereferenceInvalidObject());
  return (tria->levels[present_level]->hexes.children[present_index] != -1);
};


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
};


template <int dim>
inline
void
TriaObjectAccessor<3,dim>::operator ++ () {
  ++present_index;
				   // is index still in the range of
				   // the vector?
  while (present_index
	 >=
	 static_cast<int>(tria->levels[present_level]->hexes.hexes.size()))
    {
				       // no -> go one level up
      ++present_level;
      present_index = 0;
				       // highest level reached?
      if (present_level >= static_cast<int>(tria->levels.size()))
	{
					   // return with past the end pointer
	  present_level = present_index = -1;
	  return;
	};
    };
};


template <int dim>
inline
void
TriaObjectAccessor<3,dim>::operator -- () {
  --present_index;
				   // is index still in the range of
				   // the vector?
  while (present_index < 0) 
    {
				       // no -> go one level down
      --present_level;
				       // lowest level reached?
      if (present_level == -1) 
	{
					   // return with past the end pointer
	  present_level = present_index = -1;
	  return;
	};
				       // else
      present_index = tria->levels[present_level]->hexes.hexes.size()-1;
    };
};


/*------------------------ Functions: TriaObjectAccessor ---------------------------*/


template <int celldim, int dim>
inline
bool
TriaObjectAccessor<celldim,dim>::used () const
{
  Assert (state() == valid, ExcDereferenceInvalidObject());
  return tria->levels[present_level]->hexes.used[present_index];
};


template<int celldim, int dim>
inline
bool
TriaObjectAccessor<celldim,dim>::user_flag_set () const
{
  Assert (used(), ExcCellNotUsed());
  return tria->levels[present_level]->hexes.user_flags[present_index];
};


template<int celldim, int dim>
inline
void
TriaObjectAccessor<celldim,dim>::set_user_flag () const
{
  Assert (used(), ExcCellNotUsed());
  tria->levels[present_level]->hexes.user_flags[present_index] = true;
};


template<int celldim, int dim>
inline
void TriaObjectAccessor<celldim,dim>::clear_user_flag () const
{
  Assert (used(), ExcCellNotUsed());
  tria->levels[present_level]->hexes.user_flags[present_index] = false;
};


template<int celldim, int dim>
inline
TriaIterator<dim,TriaObjectAccessor<1,dim> >
TriaObjectAccessor<celldim,dim>::line (const unsigned int i) const
{
  Assert (used(), ExcCellNotUsed());
  Assert (i < GeometryInfo<celldim>::lines_per_cell,
	  ExcIndexRange(i,0,GeometryInfo<celldim>::lines_per_cell));

  switch(celldim)
    {
      case 2:
	    return
	      TriaIterator<dim,TriaObjectAccessor<1,dim> >
	      (
		tria,
		present_level,
		line_index (i)
	      );
      case 3:
	    if (i<4)
	      return quad(0)->line(i);
	    else
	      if (i<8)
		return quad(1)->line(i-4);
	      else
		switch (i) 
		  {
		    case 8:
			  return quad(2)->line(3);
		    case 9:
			  return quad(2)->line(1);
		    case 10:
			  return quad(4)->line(1);
		    case 11:
			  return quad(4)->line(3);
		  }
	    break;
      default:
	    Assert(false, ExcNotImplemented());
    }
  
  return TriaIterator<dim,TriaObjectAccessor<1,dim> >(tria, -1, -1, 0);
};


template<int celldim, int dim>
inline
TriaIterator<dim,TriaObjectAccessor<2,dim> >
TriaObjectAccessor<celldim,dim>::quad (const unsigned int i) const
{
  Assert (used(), ExcCellNotUsed());
  Assert (i < GeometryInfo<celldim>::quads_per_cell,
	  ExcIndexRange(i,0,GeometryInfo<celldim>::quads_per_cell));
  return
    TriaIterator<dim,TriaObjectAccessor<2,dim> >
    (
      tria,
      present_level,
      quad_index (i)
    );
};


template<int celldim, int dim>
inline
unsigned int
TriaObjectAccessor<celldim,dim>::line_index (unsigned int i) const
{
  Assert (i < GeometryInfo<celldim>::lines_per_cell,
	  ExcIndexRange(i,0,GeometryInfo<celldim>::lines_per_cell));

  switch(celldim)
    {
      case 2:
	    return tria->levels[present_level]
	      ->quads.quads[present_index].line(i);
      case 3:
	    if (i<4)
	      return quad(0)->line_index(i);
	    else
	      if (i<8)
		return quad(1)->line_index(i-4);
	      else
		switch (i) 
		  {
		    case 8:
			  return quad(2)->line_index(3);
		    case 9:
			  return quad(2)->line_index(1);
		    case 10:
			  return quad(4)->line_index(1);
		    case 11:
			  return quad(4)->line_index(3);
		  };
	    break;
      default:
	    Assert(false, ExcNotImplemented());
    }    
  return 0;
};


template<int celldim, int dim>
inline
unsigned int
TriaObjectAccessor<celldim,dim>::quad_index (unsigned int i) const
{
  Assert (i < GeometryInfo<celldim>::quads_per_cell,
	  ExcIndexRange(i,0,GeometryInfo<celldim>::quads_per_cell));

  return tria->levels[present_level]->hexes.hexes[present_index].quad(i);
};


template<int celldim, int dim>
inline
TriaIterator<dim,TriaObjectAccessor<celldim,dim> >
TriaObjectAccessor<celldim,dim>::child (const unsigned int i) const
{
  Assert (i < GeometryInfo<celldim>::children_per_cell,
	  ExcIndexRange(i,0,GeometryInfo<celldim>::children_per_cell));
  
  TriaIterator<dim,TriaObjectAccessor<celldim,dim> > q (tria, present_level+1, child_index (i));
  
#ifdef DEBUG
  if (q.state() != past_the_end)
    Assert (q->used(), ExcUnusedCellAsChild());
#endif
  return q;
};


template<int celldim, int dim>
inline
int TriaObjectAccessor<celldim,dim>::child_index (unsigned int i) const
{
  Assert (i < GeometryInfo<celldim>::children_per_cell,
	  ExcIndexRange(i,0,GeometryInfo<celldim>::children_per_cell));
  return tria->levels[present_level]->hexes.children[present_index]+i;
};


template<int celldim, int dim>
bool TriaObjectAccessor<celldim,dim>::has_children () const
{
  Assert (state() == valid, ExcDereferenceInvalidObject());
  return (tria->levels[present_level]->hexes.children[present_index] != -1);
};


template<int celldim, int dim>
inline
unsigned int
TriaObjectAccessor<celldim,dim>::max_refinement_depth () const
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
  return max (max (max (depths[0], depths[1]),
		   max (depths[2], depths[3])),
	      max (max (depths[4], depths[5]),
		   max (depths[6], depths[7])));
};


template<int celldim, int dim>
inline
void
TriaObjectAccessor<celldim,dim>::operator ++ () {
  ++present_index;
				   // is index still in the range of
				   // the vector?
  while (present_index
	 >=
	 static_cast<int>(tria->levels[present_level]->hexes.hexes.size()))
    {
				       // no -> go one level up
      ++present_level;
      present_index = 0;
				       // highest level reached?
      if (present_level >= static_cast<int>(tria->levels.size()))
	{
					   // return with past the end pointer
	  present_level = present_index = -1;
	  return;
	};
    };
};


template<int celldim, int dim>
inline
void
TriaObjectAccessor<celldim,dim>::operator -- ()
{
  --present_index;
				   // is index still in the range of
				   // the vector?
  while (present_index < 0) 
    {
				       // no -> go one level down
      --present_level;
				       // lowest level reached?
      if (present_level == -1) 
	{
					   // return with past the end pointer
	  present_level = present_index = -1;
	  return;
	};
				       // else
      present_index = tria->levels[present_level]->hexes.hexes.size()-1;
    };
};


/*------------------------ Functions: CellAccessor<dim> -----------------------*/


template <>
inline
TriaIterator<1,TriaObjectAccessor<0, 1> >
CellAccessor<1>::face (const unsigned int) const {
  Assert (false, ExcNotUsefulForThisDimension());
  return TriaIterator<1,TriaObjectAccessor<0, 1> >();
};


template <>
inline
Triangulation<2>::face_iterator
CellAccessor<2>::face (const unsigned int i) const {
  return line(i);
};


template <>
inline
Triangulation<3>::face_iterator
CellAccessor<3>::face (const unsigned int i) const {
  return quad(i);
};


template <int dim>
inline
int
CellAccessor<dim>::neighbor_index (const unsigned int i) const {
  Assert (i<GeometryInfo<dim>::faces_per_cell,
	  typename TriaAccessor<dim>::ExcInvalidNeighbor(i));
  return tria->levels[present_level]->
    neighbors[present_index*GeometryInfo<dim>::faces_per_cell+i].second;
};


template <int dim>
inline
int
CellAccessor<dim>::neighbor_level (const unsigned int i) const
{
  Assert (i<GeometryInfo<dim>::faces_per_cell,
	  typename TriaAccessor<dim>::ExcInvalidNeighbor(i));
  return tria->levels[present_level]->
    neighbors[present_index*GeometryInfo<dim>::faces_per_cell+i].first;
};



template <int dim>
inline
bool
CellAccessor<dim>::refine_flag_set () const
{
  Assert (used(), typename TriaAccessor<dim>::ExcCellNotUsed());
				   // cells flagged for refinement must be active
				   // (the @p{set_refine_flag} function checks this,
				   // but activity may change when refinement is
				   // executed and for some reason the refine
				   // flag is not cleared).
  Assert (active() ||  !tria->levels[present_level]->refine_flags[present_index],
	  ExcRefineCellNotActive());
  return tria->levels[present_level]->refine_flags[present_index];
};



template <int dim>
inline
void
CellAccessor<dim>::set_refine_flag () const
{
  Assert (used() && active(), ExcRefineCellNotActive());
  Assert (!coarsen_flag_set(),
	  ExcCellFlaggedForCoarsening());
  
  tria->levels[present_level]->refine_flags[present_index] = true;
};



template <int dim>
inline
void
CellAccessor<dim>::clear_refine_flag () const
{
  Assert (used() && active(), ExcRefineCellNotActive());
  tria->levels[present_level]->refine_flags[present_index] = false;
};



template <int dim>
inline
bool
CellAccessor<dim>::coarsen_flag_set () const
{
  Assert (used(), typename TriaAccessor<dim>::ExcCellNotUsed());
				   // cells flagged for coarsening must be active
				   // (the @p{set_refine_flag} function checks this,
				   // but activity may change when refinement is
				   // executed and for some reason the refine
				   // flag is not cleared).
  Assert (active() ||  !tria->levels[present_level]->coarsen_flags[present_index],
	  ExcRefineCellNotActive());
  return tria->levels[present_level]->coarsen_flags[present_index];
};



template <int dim>
inline
void
CellAccessor<dim>::set_coarsen_flag () const
{
  Assert (used() && active(), ExcRefineCellNotActive());
  Assert (!refine_flag_set(), ExcCellFlaggedForRefinement());
  
  tria->levels[present_level]->coarsen_flags[present_index] = true;
};



template <int dim>
inline
void
CellAccessor<dim>::clear_coarsen_flag () const
{
  Assert (used() && active(), ExcRefineCellNotActive());
  tria->levels[present_level]->coarsen_flags[present_index] = false;
};



template <int dim>
inline
TriaIterator<dim,CellAccessor<dim> >
CellAccessor<dim>::neighbor (const unsigned int i) const
{
  TriaIterator<dim,CellAccessor<dim> > q (tria, neighbor_level (i), neighbor_index (i));

#ifdef DEBUG
  if (q.state() != past_the_end)
    Assert (q->used(), typename TriaAccessor<dim>::ExcUnusedCellAsNeighbor());
#endif
  return q;
};



template <int dim>
inline
TriaIterator<dim,CellAccessor<dim> >
CellAccessor<dim>::child (const unsigned int i) const
{
  TriaIterator<dim,CellAccessor<dim> > q (tria, present_level+1, child_index (i));

#ifdef DEBUG
  if (q.state() != past_the_end)
    Assert (q->used(), typename TriaAccessor<dim>::ExcUnusedCellAsChild());
#endif
  return q;
};



template <int dim>
inline
bool
CellAccessor<dim>::active () const
{
  return !has_children();
};


#endif
