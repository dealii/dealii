/* $Id$ */
/* Copyright W. Bangerth, University of Heidelberg, 1998 */



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
distribute_local_to_global (const Vector<double> &local_source,
			    Vector<double>       &global_destination) const {
  Assert (dof_handler != 0, ExcInvalidObject());
  Assert (dof_handler->selected_fe != 0, ExcInvalidObject());
  Assert (local_source.size() == (2*dof_handler->get_fe().dofs_per_vertex +
				  dof_handler->get_fe().dofs_per_line),
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
distribute_local_to_global (const FullMatrix<double> &local_source,
			    SparseMatrix<double>     &global_destination) const {
  Assert (dof_handler != 0, ExcInvalidObject());
  Assert (dof_handler->selected_fe != 0, ExcInvalidObject());
  Assert (local_source.m() == (2*dof_handler->get_fe().dofs_per_vertex +
			       dof_handler->get_fe().dofs_per_line),
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
distribute_local_to_global (const Vector<double> &local_source,
			    Vector<double>       &global_destination) const {
  Assert (dof_handler != 0, ExcInvalidObject());
  Assert (dof_handler->selected_fe != 0, ExcInvalidObject());
  Assert (local_source.size() == (4*dof_handler->get_fe().dofs_per_vertex +
				  4*dof_handler->get_fe().dofs_per_line +
				  dof_handler->get_fe().dofs_per_quad),
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
distribute_local_to_global (const FullMatrix<double> &local_source,
			    SparseMatrix<double>     &global_destination) const {
  Assert (dof_handler != 0, ExcInvalidObject());
  Assert (dof_handler->selected_fe != 0, ExcInvalidObject());
  Assert (local_source.m() == (4*dof_handler->get_fe().dofs_per_vertex +
			       4*dof_handler->get_fe().dofs_per_line +
			       dof_handler->get_fe().dofs_per_quad),
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



/*------------------------- Functions: DoFHexAccessor -----------------------*/


template <int dim, typename BaseClass>
inline
int DoFHexAccessor<dim,BaseClass>::dof_index (const unsigned int i) const {
  Assert (dof_handler != 0, ExcInvalidObject());
				   // make sure a FE has been selected
				   // and enough room was reserved
  Assert (dof_handler->selected_fe != 0, ExcInvalidObject());
  Assert (i<dof_handler->selected_fe->dofs_per_hex,
	  ExcInvalidIndex (i, 0, dof_handler->selected_fe->dofs_per_hex));

  return dof_handler->levels[present_level]
    ->hex_dofs[present_index*dof_handler->selected_fe->dofs_per_hex+i];
};



template <int dim, typename BaseClass>
void DoFHexAccessor<dim,BaseClass>::set_dof_index (const unsigned int i,
						    const int index) const {
  Assert (dof_handler != 0, ExcInvalidObject());
				   // make sure a FE has been selected
				   // and enough room was reserved
  Assert (dof_handler->selected_fe != 0, ExcInvalidObject());
  Assert (i<dof_handler->selected_fe->dofs_per_hex,
	  ExcInvalidIndex (i, 0, dof_handler->selected_fe->dofs_per_hex));

  dof_handler->levels[present_level]
    ->hex_dofs[present_index*dof_handler->selected_fe->dofs_per_hex+i] = index;
};



template <int dim, typename BaseClass>
inline
int DoFHexAccessor<dim,BaseClass>::vertex_dof_index (const unsigned int vertex,
						      const unsigned int i) const {
  Assert (dof_handler != 0, ExcInvalidObject());
  Assert (dof_handler->selected_fe != 0, ExcInvalidObject());
  Assert (vertex<8, ExcInvalidIndex (i,0,8));
  Assert (i<dof_handler->selected_fe->dofs_per_vertex,
	  ExcInvalidIndex (i, 0, dof_handler->selected_fe->dofs_per_vertex));

  const unsigned int dof_number = (vertex_index(vertex) *
				   dof_handler->selected_fe->dofs_per_vertex +
				   i);
  return dof_handler->vertex_dofs[dof_number];
};


  
template <int dim, typename BaseClass>
void DoFHexAccessor<dim,BaseClass>::set_vertex_dof_index (const unsigned int vertex,
							   const unsigned int i,
							   const int index) const {
  Assert (dof_handler != 0, ExcInvalidObject());
  Assert (dof_handler->selected_fe != 0, ExcInvalidObject());
  Assert (vertex<8, ExcInvalidIndex (i,0,8));
  Assert (i<dof_handler->selected_fe->dofs_per_vertex,
	  ExcInvalidIndex (i, 0, dof_handler->selected_fe->dofs_per_vertex));

  const unsigned int dof_number = (vertex_index(vertex) *
				   dof_handler->selected_fe->dofs_per_vertex +
				   i);
  dof_handler->vertex_dofs[dof_number] = index;
};



template <int dim, typename BaseClass>
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
TriaIterator<dim,DoFQuadAccessor<dim,QuadAccessor<dim> > >
DoFHexAccessor<dim,BaseClass>::quad (const unsigned int i) const {
  Assert (i<6, ExcInvalidIndex (i, 0, 6));

  return TriaIterator<dim,DoFQuadAccessor<dim,QuadAccessor<dim> > >
    (
      tria,
      present_level,
      quad_index (i),
      dof_handler
    );
};



template <int dim, typename BaseClass>
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
void DoFHexAccessor<dim,BaseClass>::
distribute_local_to_global (const Vector<double> &local_source,
			    Vector<double>       &global_destination) const {
  Assert (dof_handler != 0, ExcInvalidObject());
  Assert (dof_handler->selected_fe != 0, ExcInvalidObject());
  Assert (local_source.size() == (8*dof_handler->get_fe().dofs_per_vertex +
				  12*dof_handler->get_fe().dofs_per_line +
				  6*dof_handler->get_fe().dofs_per_quad +
				  dof_handler->get_fe().dofs_per_hex),
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
void DoFHexAccessor<dim,BaseClass>::
distribute_local_to_global (const FullMatrix<double> &local_source,
			    SparseMatrix<double>     &global_destination) const {
  Assert (dof_handler != 0, ExcInvalidObject());
  Assert (dof_handler->selected_fe != 0, ExcInvalidObject());
  Assert (local_source.m() == (8*dof_handler->get_fe().dofs_per_vertex +
			       12*dof_handler->get_fe().dofs_per_line +
			       6*dof_handler->get_fe().dofs_per_quad +
			       dof_handler->get_fe().dofs_per_hex),
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
void DoFHexAccessor<dim,BaseClass>::copy_from (const DoFHexAccessor<dim,BaseClass> &a) {
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
DoFCellAccessor<1>::get_dof_values (const Vector<double> &values,
				    Vector<double>       &local_values) const {
  Assert (dof_handler != 0, ExcInvalidObject());
  Assert (&dof_handler->get_fe() != 0, ExcInvalidObject());
  Assert (local_values.size() == dof_handler->get_fe().total_dofs,
	  ExcVectorDoesNotMatch());
  Assert (values.size() == dof_handler->n_dofs(),
	  ExcVectorDoesNotMatch());
  Assert (active(), ExcNotActive());
  
  const unsigned int dofs_per_vertex = dof_handler->get_fe().dofs_per_vertex,
		     dofs_per_line   = dof_handler->get_fe().dofs_per_line;
  vector<double>::iterator next_local_value=local_values.begin();
  for (unsigned int vertex=0; vertex<2; ++vertex)
    for (unsigned int d=0; d<dofs_per_vertex; ++d)
      *next_local_value++ = values(vertex_dof_index(vertex,d));
  for (unsigned int d=0; d<dofs_per_line; ++d)
    *next_local_value++ = values(dof_index(d));
};



template <>
void
DoFCellAccessor<1>::set_dof_values (const Vector<double> &local_values,
				    Vector<double>       &values) const {
  Assert (dof_handler != 0, ExcInvalidObject());
  Assert (&dof_handler->get_fe() != 0, ExcInvalidObject());
  Assert (local_values.size() == dof_handler->get_fe().total_dofs,
	  ExcVectorDoesNotMatch());
  Assert (values.size() == dof_handler->n_dofs(),
	  ExcVectorDoesNotMatch());
  Assert (active(), ExcNotActive());
  
  const unsigned int dofs_per_vertex = dof_handler->get_fe().dofs_per_vertex,
		     dofs_per_line   = dof_handler->get_fe().dofs_per_line;
  vector<double>::const_iterator next_local_value=local_values.begin();
  for (unsigned int vertex=0; vertex<2; ++vertex)
    for (unsigned int d=0; d<dofs_per_vertex; ++d)
       values(vertex_dof_index(vertex,d)) = *next_local_value++;
  for (unsigned int d=0; d<dofs_per_line; ++d)
    values(dof_index(d)) = *next_local_value++;
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
DoFCellAccessor<2>::get_dof_values (const Vector<double> &values,
				    Vector<double>       &local_values) const {
  Assert (dof_handler != 0, ExcInvalidObject());
  Assert (&dof_handler->get_fe() != 0, ExcInvalidObject());
  Assert (local_values.size() == dof_handler->get_fe().total_dofs,
	  ExcVectorDoesNotMatch());
  Assert (values.size() == dof_handler->n_dofs(),
	  ExcVectorDoesNotMatch());
  Assert (active(), ExcNotActive());
  
  const unsigned int dofs_per_vertex = dof_handler->get_fe().dofs_per_vertex,
		     dofs_per_line   = dof_handler->get_fe().dofs_per_line,
		     dofs_per_quad   = dof_handler->get_fe().dofs_per_quad;
  vector<double>::iterator next_local_value=local_values.begin();
  for (unsigned int vertex=0; vertex<4; ++vertex)
    for (unsigned int d=0; d<dofs_per_vertex; ++d)
      *next_local_value++ = values(vertex_dof_index(vertex,d));
  for (unsigned int line=0; line<4; ++line)
    for (unsigned int d=0; d<dofs_per_line; ++d)
      *next_local_value++ = values(this->line(line)->dof_index(d));
  for (unsigned int d=0; d<dofs_per_quad; ++d)
    *next_local_value++ = values(dof_index(d));
};



template <>
void
DoFCellAccessor<2>::set_dof_values (const Vector<double> &local_values,
				    Vector<double>       &values) const {
  Assert (dof_handler != 0, ExcInvalidObject());
  Assert (&dof_handler->get_fe() != 0, ExcInvalidObject());
  Assert (local_values.size() == dof_handler->get_fe().total_dofs,
	  ExcVectorDoesNotMatch());
  Assert (values.size() == dof_handler->n_dofs(),
	  ExcVectorDoesNotMatch());
  Assert (active(), ExcNotActive());
  
  const unsigned int dofs_per_vertex = dof_handler->get_fe().dofs_per_vertex,
		     dofs_per_line   = dof_handler->get_fe().dofs_per_line,
		     dofs_per_quad   = dof_handler->get_fe().dofs_per_quad;
  vector<double>::const_iterator next_local_value=local_values.begin();
  for (unsigned int vertex=0; vertex<4; ++vertex)
    for (unsigned int d=0; d<dofs_per_vertex; ++d)
      values(vertex_dof_index(vertex,d)) = *next_local_value++;
  for (unsigned int line=0; line<4; ++line)
    for (unsigned int d=0; d<dofs_per_line; ++d)
      values(this->line(line)->dof_index(d)) = *next_local_value++;
  for (unsigned int d=0; d<dofs_per_quad; ++d)
    values(dof_index(d)) = *next_local_value++;
};



#endif



#if deal_II_dimension == 3

template <>
DoFSubstructAccessor<3>::face_iterator
DoFCellAccessor<3>::face (const unsigned int i) const {
  return quad(i);
};



template <>
void
DoFCellAccessor<3>::get_dof_values (const Vector<double> &values,
				    Vector<double>       &local_values) const {
  Assert (dof_handler != 0, ExcInvalidObject());
  Assert (&dof_handler->get_fe() != 0, ExcInvalidObject());
  Assert (local_values.size() == dof_handler->get_fe().total_dofs,
	  ExcVectorDoesNotMatch());
  Assert (values.size() == dof_handler->n_dofs(),
	  ExcVectorDoesNotMatch());
  Assert (active(), ExcNotActive());
  
  const unsigned int dofs_per_vertex = dof_handler->get_fe().dofs_per_vertex,
		     dofs_per_line   = dof_handler->get_fe().dofs_per_line,
		     dofs_per_quad   = dof_handler->get_fe().dofs_per_quad,
		     dofs_per_hex    = dof_handler->get_fe().dofs_per_hex;
  vector<double>::iterator next_local_value=local_values.begin();
  for (unsigned int vertex=0; vertex<8; ++vertex)
    for (unsigned int d=0; d<dofs_per_vertex; ++d)
      *next_local_value++ = values(vertex_dof_index(vertex,d));
  for (unsigned int line=0; line<12; ++line)
    for (unsigned int d=0; d<dofs_per_line; ++d)
      *next_local_value++ = values(this->line(line)->dof_index(d));
  for (unsigned int quad=0; quad<6; ++quad)
    for (unsigned int d=0; d<dofs_per_quad; ++d)
      *next_local_value++ = values(this->quad(quad)->dof_index(d));
  for (unsigned int d=0; d<dofs_per_hex; ++d)
    *next_local_value++ = values(dof_index(d));
};



template <>
void
DoFCellAccessor<3>::set_dof_values (const Vector<double> &local_values,
				    Vector<double>       &values) const {
  Assert (dof_handler != 0, ExcInvalidObject());
  Assert (&dof_handler->get_fe() != 0, ExcInvalidObject());
  Assert (local_values.size() == dof_handler->get_fe().total_dofs,
	  ExcVectorDoesNotMatch());
  Assert (values.size() == dof_handler->n_dofs(),
	  ExcVectorDoesNotMatch());
  Assert (active(), ExcNotActive());
  
  const unsigned int dofs_per_vertex = dof_handler->get_fe().dofs_per_vertex,
		     dofs_per_line   = dof_handler->get_fe().dofs_per_line,
		     dofs_per_quad   = dof_handler->get_fe().dofs_per_quad,
		     dofs_per_hex    = dof_handler->get_fe().dofs_per_hex;
  vector<double>::const_iterator next_local_value=local_values.begin();
  for (unsigned int vertex=0; vertex<8; ++vertex)
    for (unsigned int d=0; d<dofs_per_vertex; ++d)
      values(vertex_dof_index(vertex,d)) = *next_local_value++;
  for (unsigned int line=0; line<12; ++line)
    for (unsigned int d=0; d<dofs_per_line; ++d)
      values(this->line(line)->dof_index(d)) = *next_local_value++;
  for (unsigned int quad=0; quad<6; ++quad)
    for (unsigned int d=0; d<dofs_per_quad; ++d)
      values(this->quad(quad)->dof_index(d)) = *next_local_value++;
  for (unsigned int d=0; d<dofs_per_hex; ++d)
    values(dof_index(d)) = *next_local_value++;
};



#endif




template <int dim>
void
DoFCellAccessor<dim>::get_interpolated_dof_values (const Vector<double> &values,
						   Vector<double>       &interpolated_values) const {
  const unsigned int total_dofs = dof_handler->get_fe().total_dofs;
  
  Assert (dof_handler != 0, ExcInvalidObject());
  Assert (&dof_handler->get_fe() != 0, ExcInvalidObject());
  Assert (interpolated_values.size() == total_dofs,
	  ExcVectorDoesNotMatch());
  Assert (values.size() == dof_handler->n_dofs(),
	  ExcVectorDoesNotMatch());

  if (!has_children())
				     // if this cell has no children: simply
				     // return the exact values on this cell
    get_dof_values (values, interpolated_values);
  else
				     // otherwise clobber them from the children
    {
      Vector<double> tmp1(total_dofs);
      Vector<double> tmp2(total_dofs);
      
      interpolated_values.clear ();

      for (unsigned int child=0; child<GeometryInfo<dim>::children_per_cell;
	   ++child)
	{
					   // get the values from the present
					   // child, if necessary by
					   // interpolation itself
	  this->child(child)->get_interpolated_dof_values (values,
						     tmp1);
					   // interpolate these to the mother
					   // cell
	  dof_handler->get_fe().restrict(child).vmult (tmp2, tmp1);
	  
					   // now write those entries in tmp2
					   // which are != 0 into the output
					   // vector. Note that we may not
					   // add them up, since we would then
					   // end in adding up the contribution
					   // from nodes on boundaries of
					   // children more than once.
	  for (unsigned int i=0; i<total_dofs; ++i)
	    if (tmp2(i) != 0)
	      interpolated_values(i) = tmp2(i);
	};
    };
};



template <int dim>
void
DoFCellAccessor<dim>::set_dof_values_by_interpolation (const Vector<double> &local_values,
						       Vector<double>       &values) const {
  const unsigned int total_dofs = dof_handler->get_fe().total_dofs;
  
  Assert (dof_handler != 0, ExcInvalidObject());
  Assert (&dof_handler->get_fe() != 0, ExcInvalidObject());
  Assert (local_values.size() == total_dofs,
	  ExcVectorDoesNotMatch());
  Assert (values.size() == dof_handler->n_dofs(),
	  ExcVectorDoesNotMatch());

  if (!has_children())
				     // if this cell has no children: simply
				     // set the values on this cell
    set_dof_values (local_values, values);
  else
				     // otherwise distribute them to the children
    {
      Vector<double> tmp(total_dofs);
      
      for (unsigned int child=0; child<GeometryInfo<dim>::children_per_cell;
	   ++child)
	{
					   // prolong the given data
					   // to the present cell
	  dof_handler->get_fe().prolongate(child).vmult (tmp, local_values);
	  this->child(child)->set_dof_values_by_interpolation (tmp, values);
	};
    };
};







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
template class TriaRawIterator<2,DoFQuadAccessor<2,QuadAccessor<2> > >;
template class TriaRawIterator<2,DoFCellAccessor<2> >;
template class TriaIterator<2,DoFLineAccessor<2,LineAccessor<2> > >;
template class TriaIterator<2,DoFCellAccessor<2> >;
template class TriaActiveIterator<2,DoFLineAccessor<2,LineAccessor<2> > >;
template class TriaActiveIterator<2,DoFCellAccessor<2> >;
#endif



#if deal_II_dimension == 3
template class DoFLineAccessor<3,LineAccessor<3> >;
template class DoFQuadAccessor<3,QuadAccessor<3> >;
template class DoFHexAccessor<3,HexAccessor<3> >;
template class DoFHexAccessor<3,CellAccessor<3> >;
template class DoFCellAccessor<3>;

template class TriaRawIterator<3,DoFLineAccessor<3,LineAccessor<3> > >;
template class TriaRawIterator<3,DoFQuadAccessor<3,QuadAccessor<3> > >;
template class TriaRawIterator<3,DoFHexAccessor<3,HexAccessor<3> > >;
template class TriaRawIterator<3,DoFCellAccessor<3> >;
template class TriaIterator<3,DoFLineAccessor<3,LineAccessor<3> > >;
template class TriaIterator<3,DoFQuadAccessor<3,QuadAccessor<3> > >;
template class TriaIterator<3,DoFCellAccessor<3> >;
template class TriaActiveIterator<3,DoFLineAccessor<3,LineAccessor<3> > >;
template class TriaActiveIterator<3,DoFQuadAccessor<3,QuadAccessor<3> > >;
template class TriaActiveIterator<3,DoFCellAccessor<3> >;
#endif


