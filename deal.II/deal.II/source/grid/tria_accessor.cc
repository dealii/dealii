/* $Id$ */

#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.templates.h>

/* Note: explicit instantiations at the end of the different sections!       */



/*------------------------ Functions: LineAccessor ---------------------------*/

template <int dim>
void TriaAccessor<dim>::copy_from (const TriaAccessor<dim> &a) {
  present_level = a.present_level;
  present_index = a.present_index;
  tria = a.tria;
};

 

template <int dim>
bool TriaAccessor<dim>::operator == (const TriaAccessor<dim> &a) const {
  return ((present_level == a.present_level) &&
	  (present_index == a.present_index) &&
	  (tria == a.tria));
};



template <int dim>
bool TriaAccessor<dim>::operator != (const TriaAccessor<dim> &a) const {
  return ((present_level != a.present_level) ||
	  (present_index != a.present_index) ||
	  (tria != a.tria));
};



template <int dim>
int TriaAccessor<dim>::level () const {
  return present_level;
};


  
template <int dim>
int TriaAccessor<dim>::index () const {
  return present_index;
};
  



template <int dim>
IteratorState TriaAccessor<dim>::state () const {
  if ((present_level>=0) && (present_index>=0))
    return valid;
  else
    if ((present_level==-1) && (present_index==-1))
      return past_the_end;
    else
      return invalid;
};





template class TriaAccessor<1>;
template class TriaAccessor<2>;




/*------------------------ Functions: LineAccessor ---------------------------*/

template <int dim>
void LineAccessor<dim>::set (const Line &line) const {
  tria->levels[present_level]->lines.lines[present_index] = line;
};



template <int dim>
int LineAccessor<dim>::vertex_index (const unsigned int i) const {
  Assert (i<2, ExcInvalidIndex(i,0,1));
  return tria->levels[present_level]->lines.lines[present_index].vertex (i);
};



template <int dim>
Point<dim> &
LineAccessor<dim>::vertex (const unsigned int i) const {
  return tria->vertices[vertex_index(i)];
};



template <int dim>
bool LineAccessor<dim>::used () const {
  Assert (state() == valid, ExcDereferenceInvalidObject());
  return tria->levels[present_level]->lines.used[present_index];
};



template <int dim>
void LineAccessor<dim>::set_used_flag () const {
  Assert (state() == valid, ExcDereferenceInvalidObject());
  tria->levels[present_level]->lines.used[present_index] = true;
};
  


template <int dim>
void LineAccessor<dim>::clear_used_flag () const {
  Assert (state() == valid, ExcDereferenceInvalidObject());
  tria->levels[present_level]->lines.used[present_index] = false;
};


template <int dim>
bool LineAccessor<dim>::user_flag_set () const {
  Assert (used(), ExcRefineCellNotActive());
  return tria->levels[present_level]->lines.user_flags[present_index];
};



template <int dim>
void LineAccessor<dim>::set_user_flag () const {
  Assert (used(), ExcCellNotUsed());
  tria->levels[present_level]->lines.user_flags[present_index] = true;
};



template <int dim>
void LineAccessor<dim>::clear_user_flag () const {
  Assert (used(), ExcCellNotUsed());
  tria->levels[present_level]->lines.user_flags[present_index] = false;
};



template <int dim>
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
int LineAccessor<dim>::child_index (unsigned int i) const {
  Assert (i<2, ExcInvalidIndex(i,0,1));
  return tria->levels[present_level]->lines.children[present_index]+i;
};




template <int dim>
void LineAccessor<dim>::set_children (const int index) const {
  Assert (used(), ExcRefineCellNotUsed());
  tria->levels[present_level]->lines.children[present_index] = index;
};



template <int dim>
void LineAccessor<dim>::clear_children () const {
  set_children (-1);
};



template <int dim>
bool LineAccessor<dim>::has_children () const {
  Assert (state() == valid, ExcDereferenceInvalidObject());
  return (tria->levels[present_level]->lines.children[present_index] != -1);
}



template <int dim>
void LineAccessor<dim>::operator ++ () {
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
void LineAccessor<dim>::operator -- () {
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




template class LineAccessor<1>;
template class LineAccessor<2>;


/*------------------------ Functions: QuadAccessor ---------------------------*/

template <int dim>
void QuadAccessor<dim>::set (const Quad &quad) const {
  tria->levels[present_level]->quads.quads[present_index] = quad;
};



template <int dim>
int QuadAccessor<dim>::vertex_index (const unsigned int corner) const {
  Assert (corner<4, ExcInvalidIndex(corner,0,3));

  const int corner_convention[4] = { 0,0,1,1 };
  return line(corner)->vertex_index(corner_convention[corner]);
};



template <int dim>
Point<dim> &
QuadAccessor<dim>::vertex (const unsigned int i) const {
  return tria->vertices[vertex_index(i)];
};



template <int dim>
bool QuadAccessor<dim>::used () const {
  Assert (state() == valid, ExcDereferenceInvalidObject());
  return tria->levels[present_level]->quads.used[present_index];
};



template <int dim>
void QuadAccessor<dim>::set_used_flag () const {
  Assert (state() == valid, ExcDereferenceInvalidObject());
  tria->levels[present_level]->quads.used[present_index] = true;
};
  


template <int dim>
void QuadAccessor<dim>::clear_used_flag () const {
  Assert (state() == valid, ExcDereferenceInvalidObject());
  tria->levels[present_level]->quads.used[present_index] = false;
};




template <int dim>
bool QuadAccessor<dim>::user_flag_set () const {
  Assert (used(), ExcRefineCellNotActive());
  return tria->levels[present_level]->quads.user_flags[present_index];
};



template <int dim>
void QuadAccessor<dim>::set_user_flag () const {
  Assert (used(), ExcCellNotUsed());
  tria->levels[present_level]->quads.user_flags[present_index] = true;
};



template <int dim>
void QuadAccessor<dim>::clear_user_flag () const {
  Assert (used(), ExcCellNotUsed());
  tria->levels[present_level]->quads.user_flags[present_index] = false;
};




template <int dim>
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
unsigned int QuadAccessor<dim>::line_index (unsigned int i) const {
  Assert (i<4, ExcInvalidIndex(i,0,3));

  return tria->levels[present_level]->quads.quads[present_index].line(i);
};



template <int dim>
TriaIterator<dim,QuadAccessor<dim> >
QuadAccessor<dim>::child (const unsigned int i) const {
  Assert (i<4, ExcInvalidIndex(i,0,3));
  
  TriaIterator<2,QuadAccessor<dim> > q (tria, present_level+1, child_index (i));
  
#ifdef DEBUG
  if (q.state() != past_the_end)
    Assert (q->used(), ExcUnusedCellAsChild());
#endif
  return q;
};



template <int dim>
int QuadAccessor<dim>::child_index (unsigned int i) const {
  Assert (i<4, ExcInvalidIndex(i,0,3));
  return tria->levels[present_level]->quads.children[present_index]+i;
};



template <int dim>
void QuadAccessor<dim>::set_children (const int index) const {
  Assert (used(), ExcRefineCellNotUsed());
  tria->levels[present_level]->quads.children[present_index] = index;
};



template <int dim>
void QuadAccessor<dim>::clear_children () const {
  set_children (-1);
};



template <int dim>
bool QuadAccessor<dim>::has_children () const {
  Assert (state() == valid, ExcDereferenceInvalidObject());
  return (tria->levels[present_level]->quads.children[present_index] != -1);
};



template <int dim>
void QuadAccessor<dim>::operator ++ () {
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
void QuadAccessor<dim>::operator -- () {
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




template class QuadAccessor<2>;



/*------------------------ Functions: CellAccessor<dim> ---------------------*/

template <int dim>
int CellAccessor<dim>::neighbor_index (const unsigned int i) const {
  Assert (i<2*dim, ExcInvalidNeighbor(i));
  return tria->levels[present_level]->neighbors[present_index*2*dim+i].second;
};



template <int dim>
int CellAccessor<dim>::neighbor_level (const unsigned int i) const {
  Assert (i<2*dim, ExcInvalidNeighbor(i));
  return tria->levels[present_level]->neighbors[present_index*2*dim+i].first;
};



template <int dim>
void CellAccessor<dim>::set_neighbor (const unsigned int i,
				      const TriaIterator<dim,CellAccessor<dim> > &pointer) const {
  Assert (i<2*dim, ExcInvalidNeighbor(i));
  tria->levels[present_level]->neighbors[present_index*2*dim+i].first
    = pointer.accessor.present_level;
  tria->levels[present_level]->neighbors[present_index*2*dim+i].second
    = pointer.accessor.present_index;
};



template <int dim>
bool CellAccessor<dim>::at_boundary (const unsigned int i) const {
  Assert (used(), ExcCellNotUsed());
  Assert (i<2*dim, ExcInvalidIndex (i,0,2*dim-1));
  
  return (neighbor(i).state() != valid);
};



template <int dim>
bool CellAccessor<dim>::refine_flag_set () const {
  Assert (used() && active(), ExcRefineCellNotActive());
  return tria->levels[present_level]->refine_flags[present_index];
};



template <int dim>
void CellAccessor<dim>::set_refine_flag () const {
  Assert (used() && active(), ExcRefineCellNotActive());
  tria->levels[present_level]->refine_flags[present_index] = true;
};



template <int dim>
void CellAccessor<dim>::clear_refine_flag () const {
  Assert (used() && active(), ExcRefineCellNotActive());
  tria->levels[present_level]->refine_flags[present_index] = false;
};



template <int dim>
TriaIterator<dim,CellAccessor<dim> >
CellAccessor<dim>::neighbor (const unsigned int i) const {
  TriaIterator<dim,CellAccessor<dim> > q (tria, neighbor_level (i), neighbor_index (i));

#ifdef DEBUG
  if (q.state() != past_the_end)
    Assert (q->used(), ExcUnusedCellAsNeighbor());
#endif
  return q;
};



template <int dim>
TriaIterator<dim,CellAccessor<dim> >
CellAccessor<dim>::child (const unsigned int i) const {
  TriaIterator<dim,CellAccessor<dim> > q (tria, present_level+1, child_index (i));

#ifdef DEBUG
  if (q.state() != past_the_end)
    Assert (q->used(), ExcUnusedCellAsChild());
#endif
  return q;
};



template <int dim>
bool CellAccessor<dim>::active () const {
  return !has_children();
};



/*------------------------ Functions: CellAccessor<1> -----------------------*/


bool CellAccessor<1>::at_boundary () const {
  return at_boundary(0) || at_boundary(1);
};




/*------------------------ Functions: CellAccessor<2> -----------------------*/


bool CellAccessor<2>::at_boundary () const {
  return at_boundary(0) || at_boundary(1) || at_boundary(2) || at_boundary(3);
};











// explicit instantiations
template class CellAccessor<1>;
template class CellAccessor<2>;

template class TriaRawIterator<1,LineAccessor<1> >;
template class TriaRawIterator<1,CellAccessor<1> >;
template class TriaRawIterator<2,LineAccessor<2> >;
template class TriaRawIterator<2,QuadAccessor<2> >;
template class TriaRawIterator<2,CellAccessor<2> >;

template class TriaIterator<1,LineAccessor<1> >;
template class TriaIterator<1,CellAccessor<1> >;
template class TriaIterator<2,LineAccessor<2> >;
template class TriaIterator<2,QuadAccessor<2> >;
template class TriaIterator<2,CellAccessor<2> >;

template class TriaActiveIterator<1,LineAccessor<1> >;
template class TriaActiveIterator<1,CellAccessor<1> >;
template class TriaActiveIterator<2,LineAccessor<2> >;
template class TriaActiveIterator<2,QuadAccessor<2> >;
template class TriaActiveIterator<2,CellAccessor<2> >;


