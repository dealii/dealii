/* $Id$ */
/* Copyright W. Bangerth, University of Heidelberg, 1998 */



#include <grid/dof_accessor.h>
#include <grid/dof.h>
#include <grid/dof_levels.h>
#include <grid/tria_iterator.h>
#include <grid/tria_iterator.templates.h>
#include <fe/fe.h>

#include <lac/dvector.h>
#include <lac/dfmatrix.h>
#include <lac/dsmatrix.h>

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
	  ExcInvalidIndex (i, 0, dof_handler->selected_fe->dofs_per_line));

  return dof_handler->levels[present_level]
    ->line_dofs[present_index*dof_handler->selected_fe->dofs_per_line+i];
};





template <int dim, typename BaseClass>
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




template <int dim, typename BaseClass>
inline
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


  
template <int dim, typename BaseClass>
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



template <int dim, typename BaseClass>
void
DoFLineAccessor<dim,BaseClass>::get_dof_indices (vector<int> &dof_indices) const {
  Assert (dof_handler != 0, ExcInvalidObject());
  Assert (dof_handler->selected_fe != 0, ExcInvalidObject());
  Assert (dof_indices.size() == (2*dof_handler->get_selected_fe().dofs_per_vertex +
				 dof_handler->get_selected_fe().dofs_per_line),
	  ExcVectorDoesNotMatch());

  const unsigned int dofs_per_vertex = dof_handler->get_selected_fe().dofs_per_vertex,
		     dofs_per_line   = dof_handler->get_selected_fe().dofs_per_line;
  vector<int>::iterator next = dof_indices.begin();
  for (unsigned int vertex=0; vertex<2; ++vertex)
    for (unsigned int d=0; d<dofs_per_vertex; ++d)
      *next++ = vertex_dof_index(vertex,d);
  for (unsigned int d=0; d<dofs_per_line; ++d)
    *next++ = dof_index(d);
};



template <int dim, typename BaseClass>
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
void DoFLineAccessor<dim,BaseClass>::
distribute_local_to_global (const dVector &local_source,
			    dVector       &global_destination) const {
  Assert (dof_handler != 0, ExcInvalidObject());
  Assert (dof_handler->selected_fe != 0, ExcInvalidObject());
  Assert (local_source.size() == (2*dof_handler->get_selected_fe().dofs_per_vertex +
				  dof_handler->get_selected_fe().dofs_per_line),
	  ExcVectorDoesNotMatch());
  Assert (dof_handler->n_dofs() == global_destination.size(),
	  ExcVectorDoesNotMatch());

  const unsigned int n_dofs = local_source.size();
  
				   // get indices of dofs
  vector<int> dofs (n_dofs);
  get_dof_indices (dofs);
  
				   // distribute cell vector
  for (unsigned int j=0; j<n_dofs; ++j)
    global_destination(dofs[j]) += local_source(j);
};



template <int dim, typename BaseClass>
void DoFLineAccessor<dim,BaseClass>::
distribute_local_to_global (const dFMatrix &local_source,
			    dSMatrix       &global_destination) const {
  Assert (dof_handler != 0, ExcInvalidObject());
  Assert (dof_handler->selected_fe != 0, ExcInvalidObject());
  Assert (local_source.m() == (2*dof_handler->get_selected_fe().dofs_per_vertex +
			       dof_handler->get_selected_fe().dofs_per_line),
	  ExcVectorDoesNotMatch());
  Assert (local_source.m() == local_source.n(),
	  ExcMatrixDoesNotMatch());
  Assert (dof_handler->n_dofs() == global_destination.m(),
	  ExcMatrixDoesNotMatch());
  Assert (global_destination.m() == global_destination.n(),
	  ExcMatrixDoesNotMatch());

  const unsigned int n_dofs = local_source.m();

				   // get indices of dofs
  vector<int> dofs (n_dofs);
  get_dof_indices (dofs);
  
				   // distribute cell matrix
  for (unsigned int i=0; i<n_dofs; ++i)
    for (unsigned int j=0; j<n_dofs; ++j)
      global_destination.add(dofs[i], dofs[j], local_source(i,j));
};



template <int dim, typename BaseClass>
void DoFLineAccessor<dim,BaseClass>::copy_from (const DoFLineAccessor<dim,BaseClass> &a) {
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
	  ExcInvalidIndex (i, 0, dof_handler->selected_fe->dofs_per_quad));

  return dof_handler->levels[present_level]
    ->quad_dofs[present_index*dof_handler->selected_fe->dofs_per_quad+i];
};



template <int dim, typename BaseClass>
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



template <int dim, typename BaseClass>
inline
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


  
template <int dim, typename BaseClass>
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



template <int dim, typename BaseClass>
void
DoFQuadAccessor<dim,BaseClass>::get_dof_indices (vector<int> &dof_indices) const {
  Assert (dof_handler != 0, ExcInvalidObject());
  Assert (dof_handler->selected_fe != 0, ExcInvalidObject());
  Assert (dof_indices.size() == (4*dof_handler->get_selected_fe().dofs_per_vertex +
				 4*dof_handler->get_selected_fe().dofs_per_line +
				 dof_handler->get_selected_fe().dofs_per_quad),
	  ExcVectorDoesNotMatch());

  const unsigned int dofs_per_vertex = dof_handler->get_selected_fe().dofs_per_vertex,
		     dofs_per_line   = dof_handler->get_selected_fe().dofs_per_line,
		     dofs_per_quad   = dof_handler->get_selected_fe().dofs_per_quad;
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



template <int dim, typename BaseClass>
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
void DoFQuadAccessor<dim,BaseClass>::
distribute_local_to_global (const dVector &local_source,
			    dVector       &global_destination) const {
  Assert (dof_handler != 0, ExcInvalidObject());
  Assert (dof_handler->selected_fe != 0, ExcInvalidObject());
  Assert (local_source.size() == (4*dof_handler->get_selected_fe().dofs_per_vertex +
				  4*dof_handler->get_selected_fe().dofs_per_line +
				  dof_handler->get_selected_fe().dofs_per_quad),
	  ExcVectorDoesNotMatch());
  Assert (dof_handler->n_dofs() == global_destination.size(),
	  ExcVectorDoesNotMatch());

  const unsigned int n_dofs = local_source.size();

				   // get indices of dofs
  vector<int> dofs (n_dofs);
  get_dof_indices (dofs);
  
				   // distribute cell vector
  for (unsigned int j=0; j<n_dofs; ++j)
    global_destination(dofs[j]) += local_source(j);
};



template <int dim, typename BaseClass>
void DoFQuadAccessor<dim,BaseClass>::
distribute_local_to_global (const dFMatrix &local_source,
			    dSMatrix       &global_destination) const {
  Assert (dof_handler != 0, ExcInvalidObject());
  Assert (dof_handler->selected_fe != 0, ExcInvalidObject());
  Assert (local_source.m() == (4*dof_handler->get_selected_fe().dofs_per_vertex +
			       4*dof_handler->get_selected_fe().dofs_per_line +
			       dof_handler->get_selected_fe().dofs_per_quad),
	  ExcMatrixDoesNotMatch());
  Assert (local_source.m() == local_source.n(),
	  ExcMatrixDoesNotMatch());
  Assert (dof_handler->n_dofs() == global_destination.m(),
	  ExcMatrixDoesNotMatch());
  Assert (global_destination.m() == global_destination.n(),
	  ExcMatrixDoesNotMatch());
  
  const unsigned int n_dofs = local_source.m();

				   // get indices of dofs
  vector<int> dofs (n_dofs);
  get_dof_indices (dofs);
  
				   // distribute cell matrix
  for (unsigned int i=0; i<n_dofs; ++i)
    for (unsigned int j=0; j<n_dofs; ++j)
      global_destination.add(dofs[i], dofs[j], local_source(i,j));
};



template <int dim, typename BaseClass>
void DoFQuadAccessor<dim,BaseClass>::copy_from (const DoFQuadAccessor<dim,BaseClass> &a) {
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



#if deal_II_dimension == 1

template <>
DoFSubstructAccessor<1>::face_iterator
DoFCellAccessor<1>::face (const unsigned int) const {
  Assert (false, ExcNotUsefulForThisDimension());
  return 0;
};



template <>
void
DoFCellAccessor<1>::get_dof_values (const dVector &values,
				    dVector       &dof_values) const {
  Assert (dof_handler != 0, ExcInvalidObject());
  Assert (&dof_handler->get_selected_fe() != 0, ExcInvalidObject());
  Assert (dof_values.size() == dof_handler->get_selected_fe().total_dofs,
	  ExcVectorDoesNotMatch());
  Assert (values.size() == dof_handler->n_dofs(),
	  ExcVectorDoesNotMatch());

  const unsigned int dofs_per_vertex = dof_handler->get_selected_fe().dofs_per_vertex,
		     dofs_per_line   = dof_handler->get_selected_fe().dofs_per_line;
  vector<double>::iterator next_dof_value=dof_values.begin();
  for (unsigned int vertex=0; vertex<2; ++vertex)
    for (unsigned int d=0; d<dofs_per_vertex; ++d)
      *next_dof_value++ = values(vertex_dof_index(vertex,d));
  for (unsigned int d=0; d<dofs_per_line; ++d)
    *next_dof_value++ = values(dof_index(d));
};

#endif



#if deal_II_dimension == 2

template <>
DoFSubstructAccessor<2>::face_iterator
DoFCellAccessor<2>::face (const unsigned int i) const {
  return line(i);
};



template <>
void
DoFCellAccessor<2>::get_dof_values (const dVector &values,
				    dVector       &dof_values) const {
  Assert (dof_handler != 0, ExcInvalidObject());
  Assert (&dof_handler->get_selected_fe() != 0, ExcInvalidObject());
  Assert (dof_values.size() == dof_handler->get_selected_fe().total_dofs,
	  ExcVectorDoesNotMatch());
  Assert (values.size() == dof_handler->n_dofs(),
	  ExcVectorDoesNotMatch());

  const unsigned int dofs_per_vertex = dof_handler->get_selected_fe().dofs_per_vertex,
		     dofs_per_line   = dof_handler->get_selected_fe().dofs_per_line,
		     dofs_per_quad   = dof_handler->get_selected_fe().dofs_per_quad;
  vector<double>::iterator next_dof_value=dof_values.begin();
  for (unsigned int vertex=0; vertex<4; ++vertex)
    for (unsigned int d=0; d<dofs_per_vertex; ++d)
      *next_dof_value++ = values(vertex_dof_index(vertex,d));
  for (unsigned int line=0; line<4; ++line)
    for (unsigned int d=0; d<dofs_per_line; ++d)
      *next_dof_value++ = values(this->line(line)->dof_index(d));
  for (unsigned int d=0; d<dofs_per_quad; ++d)
    *next_dof_value++ = values(dof_index(d));
};

#endif







// explicit instantiations
#if deal_II_dimension == 1
template class DoFLineAccessor<1,CellAccessor<1> >;
template class DoFCellAccessor<1>;

template class TriaRawIterator<1,DoFCellAccessor<1> >;
template class TriaIterator<1,DoFCellAccessor<1> >;
template class TriaActiveIterator<1,DoFCellAccessor<1> >;
#endif

#if deal_II_dimension == 2
template class DoFLineAccessor<2,LineAccessor<2> >;
template class DoFQuadAccessor<2,QuadAccessor<2> >;
template class DoFQuadAccessor<2,CellAccessor<2> >;
template class DoFCellAccessor<2>;

template class TriaRawIterator<2,DoFLineAccessor<2,LineAccessor<2> > >;
template class TriaRawIterator<2,DoFCellAccessor<2> >;
template class TriaIterator<2,DoFLineAccessor<2,LineAccessor<2> > >;
template class TriaIterator<2,DoFCellAccessor<2> >;
template class TriaActiveIterator<2,DoFLineAccessor<2,LineAccessor<2> > >;
template class TriaActiveIterator<2,DoFCellAccessor<2> >;
#endif


