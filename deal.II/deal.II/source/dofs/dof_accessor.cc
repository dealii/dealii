/* $Id$ */


#include <grid/dof_accessor.h>
#include <grid/dof.h>
#include <grid/tria_iterator.h>
#include <grid/tria_iterator.templates.h>
#include <fe/fe.h>
#include <vector>





/*------------------------- Functions: DoFLineAccessor -----------------------*/


template <int dim, class BaseClass>
int DoFLineAccessor<dim,BaseClass>::dof_index (const unsigned int i) const {
  Assert (dof_handler != 0, ExcInvalidObject());
				   // make sure a FE has been selected
				   // and enough room was reserved
  Assert (dof_handler->selected_fe != 0, ExcInvalidObject());
  Assert (i<dof_handler->selected_fe->dofs_per_line,
	  ExcInvalidIndex (i, 0, dof_handler->selected_fe->dofs_per_line));

  return dof_handler->levels[present_level]
    ->line_dofs[present_index*dof_handler->selected_fe->dofs_per_line+i];
};





template <int dim, class BaseClass>
void DoFLineAccessor<dim,BaseClass>::set_dof_index (const unsigned int i,
						    const int index) const {
  Assert (dof_handler != 0, ExcInvalidObject());
				   // make sure a FE has been selected
				   // and enough room was reserved
  Assert (dof_handler->selected_fe != 0, ExcInvalidObject());
  Assert (i<dof_handler->selected_fe->dofs_per_line,
	  ExcInvalidIndex (i, 0, dof_handler->selected_fe->dofs_per_line));

  dof_handler->levels[present_level]
    ->line_dofs[present_index*dof_handler->selected_fe->dofs_per_line+i] = index;
};




template <int dim, class BaseClass>
int DoFLineAccessor<dim,BaseClass>::vertex_dof_index (const unsigned int vertex,
						      const unsigned int i) const {
  Assert (dof_handler != 0, ExcInvalidObject());
  Assert (dof_handler->selected_fe != 0, ExcInvalidObject());
  Assert (vertex<2, ExcInvalidIndex (i,0,2));
  Assert (i<dof_handler->selected_fe->dofs_per_vertex,
	  ExcInvalidIndex (i, 0, dof_handler->selected_fe->dofs_per_vertex));

  const unsigned int dof_number = (vertex_index(vertex) *
				   dof_handler->selected_fe->dofs_per_vertex +
				   i);
  return dof_handler->vertex_dofs[dof_number];
};


  
template <int dim, class BaseClass>
void DoFLineAccessor<dim,BaseClass>::set_vertex_dof_index (const unsigned int vertex,
							   const unsigned int i,
							   const int index) const {
  Assert (dof_handler != 0, ExcInvalidObject());
  Assert (dof_handler->selected_fe != 0, ExcInvalidObject());
  Assert (vertex<2, ExcInvalidIndex (i,0,2));
  Assert (i<dof_handler->selected_fe->dofs_per_vertex,
	  ExcInvalidIndex (i, 0, dof_handler->selected_fe->dofs_per_vertex));

  const unsigned int dof_number = (vertex_index(vertex) *
				   dof_handler->selected_fe->dofs_per_vertex +
				   i);
  dof_handler->vertex_dofs[dof_number] = index;
};



template <int dim, class BaseClass>
void
DoFLineAccessor<dim,BaseClass>::get_dof_indices (vector<int> &dof_indices) const {
  Assert (dof_handler != 0, ExcInvalidObject());
  Assert (&dof_handler->get_selected_fe() != 0, ExcInvalidObject());
  Assert (dof_indices.size() == 0, ExcVectorNotEmpty());
  
  dof_indices.reserve (dof_handler->get_selected_fe().total_dofs);
  for (unsigned int vertex=0; vertex<2; ++vertex)
    for (unsigned int d=0; d<dof_handler->get_selected_fe().dofs_per_vertex; ++d)
      dof_indices.push_back (vertex_dof_index(vertex,d));
  for (unsigned int d=0; d<dof_handler->get_selected_fe().dofs_per_line; ++d)
    dof_indices.push_back (dof_index(d));
};



template <int dim, class BaseClass>
double
DoFLineAccessor<dim,BaseClass>::diameter () const {
  return sqrt((vertex(1)-vertex(0)).square());
};



template <int dim, class BaseClass>
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



template <int dim, class BaseClass>
void DoFLineAccessor<dim,BaseClass>::copy_from (const DoFLineAccessor<dim,BaseClass> &a) {
  BaseClass::copy_from (a);
  set_dof_handler (a.dof_handler);
};



/*------------------------- Functions: DoFQuadAccessor -----------------------*/


template <int dim, class BaseClass>
int DoFQuadAccessor<dim,BaseClass>::dof_index (const unsigned int i) const {
  Assert (dof_handler != 0, ExcInvalidObject());
				   // make sure a FE has been selected
				   // and enough room was reserved
  Assert (dof_handler->selected_fe != 0, ExcInvalidObject());
  Assert (i<dof_handler->selected_fe->dofs_per_quad,
	  ExcInvalidIndex (i, 0, dof_handler->selected_fe->dofs_per_quad));

  return dof_handler->levels[present_level]
    ->quad_dofs[present_index*dof_handler->selected_fe->dofs_per_quad+i];
};



template <int dim, class BaseClass>
void DoFQuadAccessor<dim,BaseClass>::set_dof_index (const unsigned int i,
						    const int index) const {
  Assert (dof_handler != 0, ExcInvalidObject());
				   // make sure a FE has been selected
				   // and enough room was reserved
  Assert (dof_handler->selected_fe != 0, ExcInvalidObject());
  Assert (i<dof_handler->selected_fe->dofs_per_quad,
	  ExcInvalidIndex (i, 0, dof_handler->selected_fe->dofs_per_quad));

  dof_handler->levels[present_level]
    ->quad_dofs[present_index*dof_handler->selected_fe->dofs_per_quad+i] = index;
};



template <int dim, class BaseClass>
int DoFQuadAccessor<dim,BaseClass>::vertex_dof_index (const unsigned int vertex,
						      const unsigned int i) const {
  Assert (dof_handler != 0, ExcInvalidObject());
  Assert (dof_handler->selected_fe != 0, ExcInvalidObject());
  Assert (vertex<4, ExcInvalidIndex (i,0,4));
  Assert (i<dof_handler->selected_fe->dofs_per_vertex,
	  ExcInvalidIndex (i, 0, dof_handler->selected_fe->dofs_per_vertex));

  const unsigned int dof_number = (vertex_index(vertex) *
				   dof_handler->selected_fe->dofs_per_vertex +
				   i);
  return dof_handler->vertex_dofs[dof_number];
};


  
template <int dim, class BaseClass>
void DoFQuadAccessor<dim,BaseClass>::set_vertex_dof_index (const unsigned int vertex,
							   const unsigned int i,
							   const int index) const {
  Assert (dof_handler != 0, ExcInvalidObject());
  Assert (dof_handler->selected_fe != 0, ExcInvalidObject());
  Assert (vertex<4, ExcInvalidIndex (i,0,4));
  Assert (i<dof_handler->selected_fe->dofs_per_vertex,
	  ExcInvalidIndex (i, 0, dof_handler->selected_fe->dofs_per_vertex));

  const unsigned int dof_number = (vertex_index(vertex) *
				   dof_handler->selected_fe->dofs_per_vertex +
				   i);
  dof_handler->vertex_dofs[dof_number] = index;
};



template <int dim, class BaseClass>
void
DoFQuadAccessor<dim,BaseClass>::get_dof_indices (vector<int> &dof_indices) const {
  Assert (dof_handler != 0, ExcInvalidObject());
  Assert (&dof_handler->get_selected_fe() != 0, ExcInvalidObject());
  Assert (dof_indices.size() == 0, ExcVectorNotEmpty());
  
  dof_indices.reserve (dof_handler->get_selected_fe().total_dofs);
  for (unsigned int vertex=0; vertex<4; ++vertex)
    for (unsigned int d=0; d<dof_handler->get_selected_fe().dofs_per_vertex; ++d)
      dof_indices.push_back (vertex_dof_index(vertex,d));
  for (unsigned int line=0; line<4; ++line)
    for (unsigned int d=0; d<dof_handler->get_selected_fe().dofs_per_line; ++d)
      dof_indices.push_back (this->line(line)->dof_index(d));
  for (unsigned int d=0; d<dof_handler->get_selected_fe().dofs_per_quad; ++d)
    dof_indices.push_back (dof_index(d));
};



template <int dim, class BaseClass>
double
DoFQuadAccessor<dim,BaseClass>::diameter () const {
  return sqrt(max((vertex(2)-vertex(0)).square(),
		  (vertex(3)-vertex(1)).square()));
};




template <int dim, class BaseClass>
TriaIterator<dim,DoFLineAccessor<dim,LineAccessor<dim> > >
DoFQuadAccessor<dim,BaseClass>::line (const unsigned int i) const {
  Assert (i<4, ExcInvalidIndex (i, 0, 4));

  return TriaIterator<dim,DoFLineAccessor<dim,LineAccessor<dim> > >
    (
      tria,
      present_level,
      line_index (i),
      dof_handler
    );
};



template <int dim, class BaseClass>
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



template <int dim, class BaseClass>
void DoFQuadAccessor<dim,BaseClass>::copy_from (const DoFQuadAccessor<dim,BaseClass> &a) {
  Assert (a.dof_handler != 0, ExcInvalidObject());

  BaseClass::copy_from (a);
  set_dof_handler (a.dof_handler);
};



/*------------------------- Functions: DoFCellAccessor -----------------------*/

template <int dim>
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




DoFSubstructAccessor<1>::face_iterator
DoFCellAccessor<1>::face (const unsigned int) const {
  Assert (false, ExcNotUsefulForThisDimension());
  return 0;
};



DoFSubstructAccessor<2>::face_iterator
DoFCellAccessor<2>::face (const unsigned int i) const {
  return line(i);
};



void
DoFCellAccessor<1>::get_dof_values (const dVector  &values,
				    vector<double> &dof_values) const {
  Assert (dof_handler != 0, ExcInvalidObject());
  Assert (&dof_handler->get_selected_fe() != 0, ExcInvalidObject());
  Assert (dof_values.size() == dof_handler->get_selected_fe().total_dofs,
	  ExcVectorDoesNotMatch());
  Assert (values.n() == dof_handler->n_dofs(),
	  ExcVectorDoesNotMatch());

  vector<double>::iterator next_dof_value=dof_values.begin();
  for (unsigned int vertex=0; vertex<2; ++vertex)
    for (unsigned int d=0; d<dof_handler->get_selected_fe().dofs_per_vertex; ++d)
      *next_dof_value++ = values(vertex_dof_index(vertex,d));
  for (unsigned int d=0; d<dof_handler->get_selected_fe().dofs_per_line; ++d)
    *next_dof_value++ = values(dof_index(d));
};




void
DoFCellAccessor<2>::get_dof_values (const dVector  &values,
				    vector<double> &dof_values) const {
  Assert (dof_handler != 0, ExcInvalidObject());
  Assert (&dof_handler->get_selected_fe() != 0, ExcInvalidObject());
  Assert (dof_values.size() == dof_handler->get_selected_fe().total_dofs,
	  ExcVectorDoesNotMatch());
  Assert (values.n() == dof_handler->n_dofs(),
	  ExcVectorDoesNotMatch());

  vector<double>::iterator next_dof_value=dof_values.begin();
  for (unsigned int vertex=0; vertex<4; ++vertex)
    for (unsigned int d=0; d<dof_handler->get_selected_fe().dofs_per_vertex; ++d)
      *next_dof_value++ = values(vertex_dof_index(vertex,d));
  for (unsigned int line=0; line<4; ++line)
    for (unsigned int d=0; d<dof_handler->get_selected_fe().dofs_per_line; ++d)
      *next_dof_value++ = values(this->line(line)->dof_index(d));
  for (unsigned int d=0; d<dof_handler->get_selected_fe().dofs_per_quad; ++d)
    *next_dof_value++ = values(dof_index(d));
};









// explicit instantiations
template class DoFLineAccessor<1,CellAccessor<1> >;
template class DoFCellAccessor<1>;

template class DoFLineAccessor<2,LineAccessor<2> >;
template class DoFQuadAccessor<2,QuadAccessor<2> >;
template class DoFQuadAccessor<2,CellAccessor<2> >;
template class DoFCellAccessor<2>;

template class TriaRawIterator<1,DoFCellAccessor<1> >;
template class TriaRawIterator<2,DoFLineAccessor<2,LineAccessor<2> > >;
template class TriaRawIterator<2,DoFCellAccessor<2> >;

template class TriaIterator<1,DoFCellAccessor<1> >;
template class TriaIterator<2,DoFLineAccessor<2,LineAccessor<2> > >;
template class TriaIterator<2,DoFCellAccessor<2> >;

template class TriaActiveIterator<1,DoFCellAccessor<1> >;
template class TriaActiveIterator<2,DoFLineAccessor<2,LineAccessor<2> > >;
template class TriaActiveIterator<2,DoFCellAccessor<2> >;

