/*----------------------------   dof_accessor.templates.h     ---------------------------*/
/*      $Id$                 */
#ifndef __dof_accessor_templates_H
#define __dof_accessor_templates_H
/*----------------------------   dof_accessor.templates.h     ---------------------------*/


#include <grid/dof_accessor.h>
#include <grid/dof.h>
#include <grid/dof_levels.h>
#include <grid/tria_iterator.h>
#include <grid/tria_iterator.templates.h>
#include <fe/fe.h>

#include <lac/vector.h>
#include <lac/fullmatrix.h>
#include <lac/sparsematrix.h>

#include <vector>


/*------------------------- Functions: DoFLineAccessor -----------------------*/

template <int dim, typename BaseClass>
inline
int DoFLineAccessor<dim,BaseClass>::dof_index (const unsigned int i) const {
  Assert (dof_handler != 0, ExcInvalidObject());
				   // make sure a FE has been selected
				   // and enough room was reserved
  Assert (dof_handler->selected_fe != 0, ExcInvalidObject());
  Assert (i<dof_handler->selected_fe->dofs_per_line,
	  ExcIndexRange (i, 0, dof_handler->selected_fe->dofs_per_line));

  return dof_handler->levels[present_level]
    ->line_dofs[present_index*dof_handler->selected_fe->dofs_per_line+i];
};




template <int dim, typename BaseClass>
inline
int DoFLineAccessor<dim,BaseClass>::vertex_dof_index (const unsigned int vertex,
						      const unsigned int i) const {
  Assert (dof_handler != 0, ExcInvalidObject());
  Assert (dof_handler->selected_fe != 0, ExcInvalidObject());
  Assert (vertex<2, ExcIndexRange (i,0,2));
  Assert (i<dof_handler->selected_fe->dofs_per_vertex,
	  ExcIndexRange (i, 0, dof_handler->selected_fe->dofs_per_vertex));

  const unsigned int dof_number = (vertex_index(vertex) *
				   dof_handler->selected_fe->dofs_per_vertex +
				   i);
  return dof_handler->vertex_dofs[dof_number];
};



template <int dim, typename BaseClass>
inline
void
DoFLineAccessor<dim,BaseClass>::get_dof_indices (vector<int> &dof_indices) const {
  Assert (dof_handler != 0, ExcInvalidObject());
  Assert (dof_handler->selected_fe != 0, ExcInvalidObject());
  Assert (dof_indices.size() == (2*dof_handler->get_fe().dofs_per_vertex +
				 dof_handler->get_fe().dofs_per_line),
	  ExcVectorDoesNotMatch());

  const unsigned int dofs_per_vertex = dof_handler->get_fe().dofs_per_vertex,
		     dofs_per_line   = dof_handler->get_fe().dofs_per_line;
  vector<int>::iterator next = dof_indices.begin();
  for (unsigned int vertex=0; vertex<2; ++vertex)
    for (unsigned int d=0; d<dofs_per_vertex; ++d)
      *next++ = vertex_dof_index(vertex,d);
  for (unsigned int d=0; d<dofs_per_line; ++d)
    *next++ = dof_index(d);
};



template <int dim, typename BaseClass>
inline
TriaIterator<dim,DoFLineAccessor<dim,BaseClass> >
DoFLineAccessor<dim,BaseClass>::child (const unsigned int i) const {
  TriaIterator<dim,DoFLineAccessor<dim,BaseClass> > q (tria,
						       present_level+1,
						       child_index (i),
						       dof_handler);
  
#ifdef DEBUG
  if (q.state() != past_the_end)
    Assert (q->used(), typename TriaAccessor<dim>::ExcUnusedCellAsChild());
#endif
  return q;
};



template <int dim, typename BaseClass>
inline
void
DoFLineAccessor<dim,BaseClass>::copy_from (const DoFLineAccessor<dim,BaseClass> &a) {
  BaseClass::copy_from (a);
  set_dof_handler (a.dof_handler);
};




/*------------------------- Functions: DoFQuadAccessor -----------------------*/

template <int dim, typename BaseClass>
inline
int DoFQuadAccessor<dim,BaseClass>::dof_index (const unsigned int i) const {
  Assert (dof_handler != 0, ExcInvalidObject());
				   // make sure a FE has been selected
				   // and enough room was reserved
  Assert (dof_handler->selected_fe != 0, ExcInvalidObject());
  Assert (i<dof_handler->selected_fe->dofs_per_quad,
	  ExcIndexRange (i, 0, dof_handler->selected_fe->dofs_per_quad));

  return dof_handler->levels[present_level]
    ->quad_dofs[present_index*dof_handler->selected_fe->dofs_per_quad+i];
};



template <int dim, typename BaseClass>
inline
int DoFQuadAccessor<dim,BaseClass>::vertex_dof_index (const unsigned int vertex,
						      const unsigned int i) const {
  Assert (dof_handler != 0, ExcInvalidObject());
  Assert (dof_handler->selected_fe != 0, ExcInvalidObject());
  Assert (vertex<4, ExcIndexRange (i,0,4));
  Assert (i<dof_handler->selected_fe->dofs_per_vertex,
	  ExcIndexRange (i, 0, dof_handler->selected_fe->dofs_per_vertex));

  const unsigned int dof_number = (vertex_index(vertex) *
				   dof_handler->selected_fe->dofs_per_vertex +
				   i);
  return dof_handler->vertex_dofs[dof_number];
};


  
template <int dim, typename BaseClass>
inline
void
DoFQuadAccessor<dim,BaseClass>::get_dof_indices (vector<int> &dof_indices) const {
  Assert (dof_handler != 0, ExcInvalidObject());
  Assert (dof_handler->selected_fe != 0, ExcInvalidObject());
  Assert (dof_indices.size() == (4*dof_handler->get_fe().dofs_per_vertex +
				 4*dof_handler->get_fe().dofs_per_line +
				 dof_handler->get_fe().dofs_per_quad),
	  ExcVectorDoesNotMatch());

  const unsigned int dofs_per_vertex = dof_handler->get_fe().dofs_per_vertex,
		     dofs_per_line   = dof_handler->get_fe().dofs_per_line,
		     dofs_per_quad   = dof_handler->get_fe().dofs_per_quad;
  vector<int>::iterator next = dof_indices.begin();
  for (unsigned int vertex=0; vertex<4; ++vertex)
    for (unsigned int d=0; d<dofs_per_vertex; ++d)
      *next++ = vertex_dof_index(vertex,d);
  for (unsigned int line=0; line<4; ++line)
    for (unsigned int d=0; d<dofs_per_line; ++d)
      *next++ = this->line(line)->dof_index(d);
  for (unsigned int d=0; d<dofs_per_quad; ++d)
    *next++ = dof_index(d);
};



template <int dim, typename BaseClass>
inline
TriaIterator<dim,DoFLineAccessor<dim,LineAccessor<dim> > >
DoFQuadAccessor<dim,BaseClass>::line (const unsigned int i) const {
  Assert (i<4, ExcIndexRange (i, 0, 4));

  return TriaIterator<dim,DoFLineAccessor<dim,LineAccessor<dim> > >
    (
      tria,
      present_level,
      line_index (i),
      dof_handler
    );
};



template <int dim, typename BaseClass>
inline
TriaIterator<dim,DoFQuadAccessor<dim,BaseClass> >
DoFQuadAccessor<dim,BaseClass>::child (const unsigned int i) const {
  TriaIterator<dim,DoFQuadAccessor<dim,BaseClass> > q (tria,
						       present_level+1,
						       child_index (i),
						       dof_handler);
  
#ifdef DEBUG
  if (q.state() != past_the_end)
    Assert (q->used(), typename TriaAccessor<dim>::ExcUnusedCellAsChild());
#endif
  return q;
};



template <int dim, typename BaseClass>
inline
void
DoFQuadAccessor<dim,BaseClass>::copy_from (const DoFQuadAccessor<dim,BaseClass> &a) {
  BaseClass::copy_from (a);
  set_dof_handler (a.dof_handler);
};




/*------------------------- Functions: DoFHexAccessor -----------------------*/


template <int dim, typename BaseClass>
inline
int DoFHexAccessor<dim,BaseClass>::dof_index (const unsigned int i) const {
  Assert (dof_handler != 0, ExcInvalidObject());
				   // make sure a FE has been selected
				   // and enough room was reserved
  Assert (dof_handler->selected_fe != 0, ExcInvalidObject());
  Assert (i<dof_handler->selected_fe->dofs_per_hex,
	  ExcIndexRange (i, 0, dof_handler->selected_fe->dofs_per_hex));

  return dof_handler->levels[present_level]
    ->hex_dofs[present_index*dof_handler->selected_fe->dofs_per_hex+i];
};



template <int dim, typename BaseClass>
inline
int DoFHexAccessor<dim,BaseClass>::vertex_dof_index (const unsigned int vertex,
						      const unsigned int i) const {
  Assert (dof_handler != 0, ExcInvalidObject());
  Assert (dof_handler->selected_fe != 0, ExcInvalidObject());
  Assert (vertex<8, ExcIndexRange (i,0,8));
  Assert (i<dof_handler->selected_fe->dofs_per_vertex,
	  ExcIndexRange (i, 0, dof_handler->selected_fe->dofs_per_vertex));

  const unsigned int dof_number = (vertex_index(vertex) *
				   dof_handler->selected_fe->dofs_per_vertex +
				   i);
  return dof_handler->vertex_dofs[dof_number];
};


  
template <int dim, typename BaseClass>
inline
void
DoFHexAccessor<dim,BaseClass>::get_dof_indices (vector<int> &dof_indices) const {
  Assert (dof_handler != 0, ExcInvalidObject());
  Assert (dof_handler->selected_fe != 0, ExcInvalidObject());
  Assert (dof_indices.size() == (8*dof_handler->get_fe().dofs_per_vertex +
				 12*dof_handler->get_fe().dofs_per_line +
				 6*dof_handler->get_fe().dofs_per_quad +
				 dof_handler->get_fe().dofs_per_hex),
	  ExcVectorDoesNotMatch());

  const unsigned int dofs_per_vertex = dof_handler->get_fe().dofs_per_vertex,
		     dofs_per_line   = dof_handler->get_fe().dofs_per_line,
		     dofs_per_quad   = dof_handler->get_fe().dofs_per_quad,
		     dofs_per_hex    = dof_handler->get_fe().dofs_per_hex;
  vector<int>::iterator next = dof_indices.begin();
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
};



template <int dim, typename BaseClass>
inline
TriaIterator<dim,DoFLineAccessor<dim,LineAccessor<dim> > >
DoFHexAccessor<dim,BaseClass>::line (const unsigned int i) const {
  TriaIterator<dim,LineAccessor<dim> > l = BaseClass::line(i);
  return TriaIterator<dim,DoFLineAccessor<dim,LineAccessor<dim> > >
    (
      tria,
      present_level,
      l->index(),
      dof_handler
    );
};



template <int dim, typename BaseClass>
inline
TriaIterator<dim,DoFQuadAccessor<dim,QuadAccessor<dim> > >
DoFHexAccessor<dim,BaseClass>::quad (const unsigned int i) const {
  Assert (i<6, ExcIndexRange (i, 0, 6));

  return TriaIterator<dim,DoFQuadAccessor<dim,QuadAccessor<dim> > >
    (
      tria,
      present_level,
      quad_index (i),
      dof_handler
    );
};



template <int dim, typename BaseClass>
inline
TriaIterator<dim,DoFHexAccessor<dim,BaseClass> >
DoFHexAccessor<dim,BaseClass>::child (const unsigned int i) const {
  TriaIterator<dim,DoFHexAccessor<dim,BaseClass> > q (tria,
						      present_level+1,
						      child_index (i),
						      dof_handler);
  
#ifdef DEBUG
  if (q.state() != past_the_end)
    Assert (q->used(), typename TriaAccessor<dim>::ExcUnusedCellAsChild());
#endif
  return q;
};



template <int dim, typename BaseClass>
void DoFHexAccessor<dim,BaseClass>::copy_from (const DoFHexAccessor<dim,BaseClass> &a) {
  BaseClass::copy_from (a);
  set_dof_handler (a.dof_handler);
};




/*------------------------- Functions: DoFCellAccessor -----------------------*/


template <int dim>
inline
TriaIterator<dim,DoFCellAccessor<dim> >
DoFCellAccessor<dim>::neighbor (const unsigned int i) const {
  TriaIterator<dim,DoFCellAccessor<dim> > q (tria,
					     neighbor_level (i),
					     neighbor_index (i),
					     dof_handler);
  
#ifdef DEBUG
  if (q.state() != past_the_end)
    Assert (q->used(), typename TriaAccessor<dim>::ExcUnusedCellAsNeighbor());
#endif
  return q;
};



template <int dim>
inline
TriaIterator<dim,DoFCellAccessor<dim> >
DoFCellAccessor<dim>::child (const unsigned int i) const {
  TriaIterator<dim,DoFCellAccessor<dim> > q (tria,
					     present_level+1,
					     child_index (i),
					     dof_handler);
  
#ifdef DEBUG
  if (q.state() != past_the_end)
    Assert (q->used(), typename TriaAccessor<dim>::ExcUnusedCellAsChild());
#endif
  return q;
};







/*----------------------------   dof_accessor.templates.h     ---------------------------*/
/* end of #ifndef __dof_accessor_templates_H */
#endif
/*----------------------------   dof_accessor.templates.h     ---------------------------*/
