/*----------------------------   tria_accessor.templates.h     ---------------------------*/
/*      $Id$                 */
#ifndef __tria_accessor_templates_H
#define __tria_accessor_templates_H
/*----------------------------   tria_accessor.templates.h     ---------------------------*/


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
LineAccessor<dim>::used () const {
  Assert (state() == valid, ExcDereferenceInvalidObject());
  return tria->levels[present_level]->lines.used[present_index];
};




template <int dim>
inline
bool
LineAccessor<dim>::user_flag_set () const {
  Assert (used(), ExcCellNotUsed());
  return tria->levels[present_level]->lines.user_flags[present_index];
};



template <int dim>
inline
void
LineAccessor<dim>::set_user_flag () const {
  Assert (used(), ExcCellNotUsed());
  tria->levels[present_level]->lines.user_flags[present_index] = true;
};



template <int dim>
inline
void
LineAccessor<dim>::clear_user_flag () const {
  Assert (used(), ExcCellNotUsed());
  tria->levels[present_level]->lines.user_flags[present_index] = false;
};




template <int dim>
inline
TriaIterator<dim,LineAccessor<dim> >
LineAccessor<dim>::child (const unsigned int i) const {
  Assert (i<2, ExcInvalidIndex(i,0,1));
  
  TriaIterator<dim,LineAccessor<dim> > q (tria, present_level+1, child_index (i));
  
#ifdef DEBUG
  if (q.state() != past_the_end)
    Assert (q->used(), ExcUnusedCellAsChild());
#endif
  return q;
};



template <int dim>
inline
int
LineAccessor<dim>::child_index (unsigned int i) const {
  Assert (i<2, ExcInvalidIndex(i,0,1));
  return tria->levels[present_level]->lines.children[present_index]+i;
};




template <int dim>
inline
bool
LineAccessor<dim>::has_children () const {
  Assert (state() == valid, ExcDereferenceInvalidObject());
  return (tria->levels[present_level]->lines.children[present_index] != -1);
}



template <int dim>
inline
unsigned int
LineAccessor<dim>::max_refinement_depth () const
{
  if (!has_children())
    return 0;

  const unsigned int depths[2] = { child(0)->max_refinement_depth() + 1,
				   child(1)->max_refinement_depth() + 1  };
  return max (depths[0], depths[1]);
};
      
	    


template <int dim>
inline
void
LineAccessor<dim>::operator ++ () {
  ++present_index;
				   // is index still in the range of
				   // the vector?
  while (present_index
	 >=
	 (int)tria->levels[present_level]->lines.lines.size())
    {
				       // no -> go one level up until we find
				       // one with more than zero cells
      ++present_level;
      present_index = 0;
				       // highest level reached?
      if (present_level >= (int)tria->levels.size()) 
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
LineAccessor<dim>::operator -- () {
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
QuadAccessor<dim>::used () const {
  Assert (state() == valid, ExcDereferenceInvalidObject());
  return tria->levels[present_level]->quads.used[present_index];
};




template <int dim>
inline
bool
QuadAccessor<dim>::user_flag_set () const {
  Assert (used(), ExcCellNotUsed());
  return tria->levels[present_level]->quads.user_flags[present_index];
};



template <int dim>
inline
void
QuadAccessor<dim>::set_user_flag () const {
  Assert (used(), ExcCellNotUsed());
  tria->levels[present_level]->quads.user_flags[present_index] = true;
};



template <int dim>
inline
void
QuadAccessor<dim>::clear_user_flag () const {
  Assert (used(), ExcCellNotUsed());
  tria->levels[present_level]->quads.user_flags[present_index] = false;
};




template <int dim>
inline
TriaIterator<dim,LineAccessor<dim> >
QuadAccessor<dim>::line (const unsigned int i) const {
  return
    TriaIterator<dim,LineAccessor<dim> >
    (
      tria,
      present_level,
      line_index (i)
    );
};



template <int dim>
inline
unsigned int
QuadAccessor<dim>::line_index (unsigned int i) const {
  Assert (i<4, ExcInvalidIndex(i,0,3));

  return tria->levels[present_level]->quads.quads[present_index].line(i);
};



template <int dim>
inline
TriaIterator<dim,QuadAccessor<dim> >
QuadAccessor<dim>::child (const unsigned int i) const {
  Assert (i<4, ExcInvalidIndex(i,0,3));
  
  TriaIterator<dim,QuadAccessor<dim> > q (tria, present_level+1, child_index (i));
  
#ifdef DEBUG
  if (q.state() != past_the_end)
    Assert (q->used(), ExcUnusedCellAsChild());
#endif
  return q;
};



template <int dim>
inline
int QuadAccessor<dim>::child_index (unsigned int i) const {
  Assert (i<4, ExcInvalidIndex(i,0,3));
  return tria->levels[present_level]->quads.children[present_index]+i;
};




template <int dim>
inline
bool
QuadAccessor<dim>::has_children () const {
  Assert (state() == valid, ExcDereferenceInvalidObject());
  return (tria->levels[present_level]->quads.children[present_index] != -1);
};



template <int dim>
inline
unsigned int
QuadAccessor<dim>::max_refinement_depth () const
{
  if (!has_children())
    return 0;

  const unsigned int depths[4] = { child(0)->max_refinement_depth() + 1,
				   child(1)->max_refinement_depth() + 1,
				   child(2)->max_refinement_depth() + 1,
				   child(3)->max_refinement_depth() + 1 };
  return max (max (depths[0], depths[1]),
	      max (depths[2], depths[3]));
};



template <int dim>
inline
void
QuadAccessor<dim>::operator ++ () {
  ++present_index;
				   // is index still in the range of
				   // the vector?
  while (present_index
	 >=
	 (int)tria->levels[present_level]->quads.quads.size()) 
    {
				       // no -> go one level up
      ++present_level;
      present_index = 0;
				       // highest level reached?
      if (present_level >= (int)tria->levels.size()) 
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
QuadAccessor<dim>::operator -- () {
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
HexAccessor<dim>::used () const {
  Assert (state() == valid, ExcDereferenceInvalidObject());
  return tria->levels[present_level]->hexes.used[present_index];
};



template <int dim>
inline
bool
HexAccessor<dim>::user_flag_set () const {
  Assert (used(), ExcCellNotUsed());
  return tria->levels[present_level]->hexes.user_flags[present_index];
};



template <int dim>
inline
void
HexAccessor<dim>::set_user_flag () const {
  Assert (used(), ExcCellNotUsed());
  tria->levels[present_level]->hexes.user_flags[present_index] = true;
};



template <int dim>
inline
void HexAccessor<dim>::clear_user_flag () const {
  Assert (used(), ExcCellNotUsed());
  tria->levels[present_level]->hexes.user_flags[present_index] = false;
};




template <int dim>
inline
TriaIterator<dim,LineAccessor<dim> >
HexAccessor<dim>::line (const unsigned int i) const {
  Assert (used(), ExcCellNotUsed());
  Assert (i<12, ExcInvalidIndex (i,0,11));
  
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
	};
  Assert (false, ExcInvalidIndex(i,0,11));
  return TriaIterator<dim,LineAccessor<dim> >(tria, -1, -1, 0);
};



template <int dim>
inline
TriaIterator<dim,QuadAccessor<dim> >
HexAccessor<dim>::quad (const unsigned int i) const {
  Assert (used(), ExcCellNotUsed());
  return
    TriaIterator<dim,QuadAccessor<dim> >
    (
      tria,
      present_level,
      quad_index (i)
    );
};



template <int dim>
inline
unsigned int
HexAccessor<dim>::line_index (unsigned int i) const {
  Assert (i<12, ExcInvalidIndex(i,0,11));

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
  Assert (false, ExcInvalidIndex(i,0,11));
  return 0;
};



template <int dim>
inline
unsigned int
HexAccessor<dim>::quad_index (unsigned int i) const {
  Assert (i<6, ExcInvalidIndex(i,0,5));

  return tria->levels[present_level]->hexes.hexes[present_index].quad(i);
};



template <int dim>
inline
TriaIterator<dim,HexAccessor<dim> >
HexAccessor<dim>::child (const unsigned int i) const {
  Assert (i<6, ExcInvalidIndex(i,0,5));
  
  TriaIterator<dim,HexAccessor<dim> > q (tria, present_level+1, child_index (i));
  
#ifdef DEBUG
  if (q.state() != past_the_end)
    Assert (q->used(), ExcUnusedCellAsChild());
#endif
  return q;
};



template <int dim>
inline
int HexAccessor<dim>::child_index (unsigned int i) const {
  Assert (i<8, ExcInvalidIndex(i,0,7));
  return tria->levels[present_level]->hexes.children[present_index]+i;
};



template <int dim>
bool HexAccessor<dim>::has_children () const {
  Assert (state() == valid, ExcDereferenceInvalidObject());
  return (tria->levels[present_level]->hexes.children[present_index] != -1);
};



template <int dim>
inline
unsigned int
HexAccessor<dim>::max_refinement_depth () const
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



template <int dim>
inline
void
HexAccessor<dim>::operator ++ () {
  ++present_index;
				   // is index still in the range of
				   // the vector?
  while (present_index
	 >=
	 (int)tria->levels[present_level]->hexes.hexes.size()) 
    {
				       // no -> go one level up
      ++present_level;
      present_index = 0;
				       // highest level reached?
      if (present_level >= (int)tria->levels.size()) 
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
HexAccessor<dim>::operator -- () {
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
Triangulation<1>::face_iterator
CellAccessor<1>::face (const unsigned int) const {
  Assert (false, ExcNotUsefulForThisDimension());
  return 0;
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
	  typename TriaSubstructAccessor<dim>::ExcInvalidNeighbor(i));
  return tria->levels[present_level]->
    neighbors[present_index*GeometryInfo<dim>::faces_per_cell+i].second;
};



template <int dim>
inline
int
CellAccessor<dim>::neighbor_level (const unsigned int i) const {
  Assert (i<GeometryInfo<dim>::faces_per_cell,
	  typename TriaSubstructAccessor<dim>::ExcInvalidNeighbor(i));
  return tria->levels[present_level]->
    neighbors[present_index*GeometryInfo<dim>::faces_per_cell+i].first;
};




template <int dim>
inline
bool
CellAccessor<dim>::refine_flag_set () const {
  Assert (used(), typename TriaSubstructAccessor<dim>::ExcCellNotUsed());
				   // cells flagged for refinement must be active
				   // (the #set_refine_flag# function checks this,
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
CellAccessor<dim>::set_refine_flag () const {
  Assert (used() && active(), ExcRefineCellNotActive());
  Assert (!coarsen_flag_set(), ExcCellFlaggedForCoarsening());
  
  tria->levels[present_level]->refine_flags[present_index] = true;
};



template <int dim>
inline
void
CellAccessor<dim>::clear_refine_flag () const {
  Assert (used() && active(), ExcRefineCellNotActive());
  tria->levels[present_level]->refine_flags[present_index] = false;
};



template <int dim>
inline
bool
CellAccessor<dim>::coarsen_flag_set () const {
  Assert (used(), typename TriaSubstructAccessor<dim>::ExcCellNotUsed());
				   // cells flagged for coarsening must be active
				   // (the #set_refine_flag# function checks this,
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
CellAccessor<dim>::set_coarsen_flag () const {
  Assert (used() && active(), ExcRefineCellNotActive());
  Assert (!refine_flag_set(), ExcCellFlaggedForRefinement());
  
  tria->levels[present_level]->coarsen_flags[present_index] = true;
};



template <int dim>
inline
void
CellAccessor<dim>::clear_coarsen_flag () const {
  Assert (used() && active(), ExcRefineCellNotActive());
  tria->levels[present_level]->coarsen_flags[present_index] = false;
};



template <int dim>
inline
TriaIterator<dim,CellAccessor<dim> >
CellAccessor<dim>::neighbor (const unsigned int i) const {
  TriaIterator<dim,CellAccessor<dim> > q (tria, neighbor_level (i), neighbor_index (i));

#ifdef DEBUG
  if (q.state() != past_the_end)
    Assert (q->used(),
	    typename TriaSubstructAccessor<dim>::ExcUnusedCellAsNeighbor());
#endif
  return q;
};



template <int dim>
inline
TriaIterator<dim,CellAccessor<dim> >
CellAccessor<dim>::child (const unsigned int i) const {
  TriaIterator<dim,CellAccessor<dim> > q (tria, present_level+1, child_index (i));

#ifdef DEBUG
  if (q.state() != past_the_end)
    Assert (q->used(),
	    typename TriaSubstructAccessor<dim>::ExcUnusedCellAsChild());
#endif
  return q;
};



template <int dim>
inline
bool
CellAccessor<dim>::active () const {
  return !has_children();
};






/*----------------------------   tria_accessor.templates.h     ---------------------------*/
/* end of #ifndef __tria_accessor_templates_H */
#endif
/*----------------------------   tria_accessor.templates.h     ---------------------------*/
