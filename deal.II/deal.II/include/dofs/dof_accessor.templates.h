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

template <int dim>
inline
int DoFObjectAccessor<1, dim>::dof_index (const unsigned int i) const {
  Assert (dof_handler != 0, ExcInvalidObject());
				   // make sure a FE has been selected
				   // and enough room was reserved
  Assert (dof_handler->selected_fe != 0, ExcInvalidObject());
  Assert (i<dof_handler->selected_fe->dofs_per_line,
	  ExcIndexRange (i, 0, dof_handler->selected_fe->dofs_per_line));

  return dof_handler->levels[present_level]
    ->line_dofs[present_index*dof_handler->selected_fe->dofs_per_line+i];
};




template <int dim>
inline
int DoFObjectAccessor<1, dim>::vertex_dof_index (const unsigned int vertex,
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



template <int dim>
inline
void
DoFObjectAccessor<1, dim>::get_dof_indices (vector<int> &dof_indices) const {
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



template <int dim>
inline
TriaIterator<dim,DoFObjectAccessor<1, dim> >
DoFObjectAccessor<1, dim>::child (const unsigned int i) const {
  TriaIterator<dim,DoFObjectAccessor<1, dim> > q (tria,
						       present_level+1,
						       child_index (i),
						       dof_handler);
  
#ifdef DEBUG
  if (q.state() != past_the_end)
    Assert (q->used(), typename TriaAccessor<dim>::ExcUnusedCellAsChild());
#endif
  return q;
};



template <int dim>
inline
void
DoFObjectAccessor<1, dim>::copy_from (const DoFObjectAccessor<1, dim> &a) {
  DoFObjectAccessor_Inheritance<1,dim>::BaseClass::copy_from (a);
  set_dof_handler (a.dof_handler);
};




/*------------------------- Functions: DoFObjectAccessor<2, dim> -----------------------*/

template <int dim>
inline
int DoFObjectAccessor<2, dim>::dof_index (const unsigned int i) const {
  Assert (dof_handler != 0, ExcInvalidObject());
				   // make sure a FE has been selected
				   // and enough room was reserved
  Assert (dof_handler->selected_fe != 0, ExcInvalidObject());
  Assert (i<dof_handler->selected_fe->dofs_per_quad,
	  ExcIndexRange (i, 0, dof_handler->selected_fe->dofs_per_quad));

  return dof_handler->levels[present_level]
    ->quad_dofs[present_index*dof_handler->selected_fe->dofs_per_quad+i];
};



template <int dim>
inline
int DoFObjectAccessor<2, dim>::vertex_dof_index (const unsigned int vertex,
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


  
template <int dim>
inline
void
DoFObjectAccessor<2, dim>::get_dof_indices (vector<int> &dof_indices) const {
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



template <int dim>
inline
TriaIterator<dim,DoFObjectAccessor<1, dim> >
DoFObjectAccessor<2, dim>::line (const unsigned int i) const {
  Assert (i<4, ExcIndexRange (i, 0, 4));

  return TriaIterator<dim,DoFObjectAccessor<1, dim> >
    (
      tria,
      present_level,
      line_index (i),
      dof_handler
    );
};



template <int dim>
inline
TriaIterator<dim,DoFObjectAccessor<2, dim> >
DoFObjectAccessor<2, dim>::child (const unsigned int i) const {
  TriaIterator<dim,DoFObjectAccessor<2, dim> > q (tria,
						  present_level+1,
						  child_index (i),
						  dof_handler);
  
#ifdef DEBUG
  if (q.state() != past_the_end)
    Assert (q->used(), typename TriaAccessor<dim>::ExcUnusedCellAsChild());
#endif
  return q;
};



template <int dim>
inline
void
DoFObjectAccessor<2, dim>::copy_from (const DoFObjectAccessor<2, dim> &a) {
  DoFObjectAccessor_Inheritance<2,dim>::BaseClass::copy_from (a);
  set_dof_handler (a.dof_handler);
};




/*------------------------- Functions: DoFObjectAccessor -----------------------*/


template <int dim>
inline
int DoFObjectAccessor<3, dim>::dof_index (const unsigned int i) const {
  Assert (dof_handler != 0, ExcInvalidObject());
				   // make sure a FE has been selected
				   // and enough room was reserved
  Assert (dof_handler->selected_fe != 0, ExcInvalidObject());
  Assert (i<dof_handler->selected_fe->dofs_per_hex,
	  ExcIndexRange (i, 0, dof_handler->selected_fe->dofs_per_hex));

  return dof_handler->levels[present_level]
    ->hex_dofs[present_index*dof_handler->selected_fe->dofs_per_hex+i];
};



template <int dim>
inline
int DoFObjectAccessor<3, dim>::vertex_dof_index (const unsigned int vertex,
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


  
template <int dim>
inline
void
DoFObjectAccessor<3, dim>::get_dof_indices (vector<int> &dof_indices) const {
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



template <int dim>
inline
TriaIterator<dim,DoFObjectAccessor<1, dim> >
DoFObjectAccessor<3, dim>::line (const unsigned int i) const {
  TriaIterator<dim,TriaObjectAccessor<1, dim> > l = DoFObjectAccessor_Inheritance<3,dim>::BaseClass::line(i);
  return TriaIterator<dim,DoFObjectAccessor<1, dim> >
    (
      tria,
      present_level,
      l->index(),
      dof_handler
    );
};



template <int dim>
inline
TriaIterator<dim,DoFObjectAccessor<2, dim> >
DoFObjectAccessor<3, dim>::quad (const unsigned int i) const {
  Assert (i<6, ExcIndexRange (i, 0, 6));

  return TriaIterator<dim,DoFObjectAccessor<2, dim> >
    (
      tria,
      present_level,
      quad_index (i),
      dof_handler
    );
};



template <int dim>
inline
TriaIterator<dim,DoFObjectAccessor<3, dim> >
DoFObjectAccessor<3, dim>::child (const unsigned int i) const {
  TriaIterator<dim,DoFObjectAccessor<3, dim> > q (tria,
						  present_level+1,
						  child_index (i),
						  dof_handler);
  
#ifdef DEBUG
  if (q.state() != past_the_end)
    Assert (q->used(), typename TriaAccessor<dim>::ExcUnusedCellAsChild());
#endif
  return q;
};



template <int dim>
void DoFObjectAccessor<3, dim>::copy_from (const DoFObjectAccessor<3, dim> &a) {
  DoFObjectAccessor_Inheritance<3,dim>::BaseClass::copy_from (a);
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
