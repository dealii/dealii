/* $Id$ */
/* Copyright W. Bangerth, University of Heidelberg, 1998 */



#include <grid/mg_dof_accessor.h>
#include <grid/mg_dof.h>
#include <grid/dof_levels.h>
#include <grid/tria_iterator.h>
#include <grid/tria_iterator.templates.h>
#include <fe/fe.h>

#include <lac/vector.h>
#include <lac/fullmatrix.h>
#include <lac/sparsematrix.h>






/* ------------------------ MGDoFLineAccessor --------------------------- */

template <int dim, typename BaseClass>
MGDoFLineAccessor<dim,BaseClass>::MGDoFLineAccessor (Triangulation<dim> *tria,
						     const int           level,
						     const int           index,
						     const AccessorData *local_data) :
		MGDoFAccessor<dim> (local_data),
		DoFLineAccessor<dim,BaseClass>(tria,level,index,local_data) {};



template <int dim, typename BaseClass>
int MGDoFLineAccessor<dim,BaseClass>::mg_dof_index (const unsigned int i) const {
  Assert (dof_handler != 0, ExcInvalidObject());
  Assert (mg_dof_handler != 0, ExcInvalidObject());
				   // make sure a FE has been selected
				   // and enough room was reserved
  Assert (&dof_handler->get_fe() != 0, ExcInvalidObject());
  Assert (i<dof_handler->get_fe().dofs_per_line,
	  ExcIndexRange (i, 0, dof_handler->get_fe().dofs_per_line));

  return mg_dof_handler->mg_levels[present_level]
    ->line_dofs[present_index*dof_handler->get_fe().dofs_per_line+i];
};





template <int dim, typename BaseClass>
void MGDoFLineAccessor<dim,BaseClass>::set_mg_dof_index (const unsigned int i,
							 const int index) const {
  Assert (dof_handler != 0, ExcInvalidObject());
  Assert (mg_dof_handler != 0, ExcInvalidObject());
				   // make sure a FE has been selected
				   // and enough room was reserved
  Assert (&dof_handler->get_fe() != 0, ExcInvalidObject());
  Assert (i<dof_handler->get_fe().dofs_per_line,
	  ExcIndexRange (i, 0, dof_handler->get_fe().dofs_per_line));

  mg_dof_handler->mg_levels[present_level]
    ->line_dofs[present_index*dof_handler->get_fe().dofs_per_line+i] = index;
};




template <int dim, typename BaseClass>
int MGDoFLineAccessor<dim,BaseClass>::mg_vertex_dof_index (const unsigned int vertex,
							   const unsigned int i) const {
  Assert (dof_handler != 0, ExcInvalidObject());
  Assert (mg_dof_handler != 0, ExcInvalidObject());
  Assert (&dof_handler->get_fe() != 0, ExcInvalidObject());
  Assert (vertex<2, ExcIndexRange (i,0,2));
  Assert (i<dof_handler->get_fe().dofs_per_vertex,
	  ExcIndexRange (i, 0, dof_handler->get_fe().dofs_per_vertex));

  return (mg_dof_handler->mg_vertex_dofs[vertex_index(vertex)]
	  .get_index (present_level, i, dof_handler->get_fe().dofs_per_vertex));
};


  
template <int dim, typename BaseClass>
void MGDoFLineAccessor<dim,BaseClass>::set_mg_vertex_dof_index (const unsigned int vertex,
								const unsigned int i,
								const int index) const {
  Assert (dof_handler != 0, ExcInvalidObject());
  Assert (mg_dof_handler != 0, ExcInvalidObject());
  Assert (&dof_handler->get_fe() != 0, ExcInvalidObject());
  Assert (vertex<2, ExcIndexRange (i,0,2));
  Assert (i<dof_handler->get_fe().dofs_per_vertex,
	  ExcIndexRange (i, 0, dof_handler->get_fe().dofs_per_vertex));

  mg_dof_handler->mg_vertex_dofs[vertex_index(vertex)]
    .set_index (present_level, i, dof_handler->get_fe().dofs_per_vertex, index);
};



template <int dim, typename BaseClass>
void
MGDoFLineAccessor<dim,BaseClass>::get_mg_dof_indices (vector<int> &dof_indices) const {
  Assert (dof_handler != 0, ExcInvalidObject());
  Assert (mg_dof_handler != 0, ExcInvalidObject());
  Assert (&dof_handler->get_fe() != 0, ExcInvalidObject());
  Assert (dof_indices.size() == (2*dof_handler->get_fe().dofs_per_vertex +
				 dof_handler->get_fe().dofs_per_line),
	  ExcVectorDoesNotMatch());

  const unsigned int dofs_per_vertex = dof_handler->get_fe().dofs_per_vertex,
		     dofs_per_line   = dof_handler->get_fe().dofs_per_line;
  vector<int>::iterator next = dof_indices.begin();
  for (unsigned int vertex=0; vertex<2; ++vertex)
    for (unsigned int d=0; d<dofs_per_vertex; ++d)
      *next++ = mg_vertex_dof_index(vertex,d);
  for (unsigned int d=0; d<dofs_per_line; ++d)
    *next++ = mg_dof_index(d);
};



template <int dim, typename BaseClass>
TriaIterator<dim,MGDoFLineAccessor<dim,BaseClass> >
MGDoFLineAccessor<dim,BaseClass>::child (const unsigned int i) const {
  TriaIterator<dim,MGDoFLineAccessor<dim,BaseClass> > q (tria,
							 present_level+1,
							 child_index (i),
							 mg_dof_handler);
  
#ifdef DEBUG
  if (q.state() != past_the_end)
    Assert (q->used(), typename TriaAccessor<dim>::ExcUnusedCellAsChild());
#endif
  return q;
};



template <int dim, typename BaseClass>
void
MGDoFLineAccessor<dim,BaseClass>::copy_from (const MGDoFLineAccessor<dim,BaseClass> &a) {
  DoFLineAccessor<dim,BaseClass>::copy_from (a);
  set_mg_dof_handler (a.mg_dof_handler);
};




/* ------------------------ MGDoFQuadAccessor --------------------------- */

template <int dim, typename BaseClass>
MGDoFQuadAccessor<dim,BaseClass>::MGDoFQuadAccessor (Triangulation<dim> *tria,
						     const int           level,
						     const int           index,
						     const AccessorData *local_data) :
		MGDoFAccessor<dim> (local_data),
		DoFQuadAccessor<dim,BaseClass>(tria,level,index,local_data) {};



template <int dim, typename BaseClass>
inline
int MGDoFQuadAccessor<dim,BaseClass>::mg_dof_index (const unsigned int i) const {
  Assert (dof_handler != 0, ExcInvalidObject());
  Assert (mg_dof_handler != 0, ExcInvalidObject());
				   // make sure a FE has been selected
				   // and enough room was reserved
  Assert (&dof_handler->get_fe() != 0, ExcInvalidObject());
  Assert (i<dof_handler->get_fe().dofs_per_quad,
	  ExcIndexRange (i, 0, dof_handler->get_fe().dofs_per_quad));

  return mg_dof_handler->mg_levels[present_level]
    ->quad_dofs[present_index*dof_handler->get_fe().dofs_per_quad+i];
};



template <int dim, typename BaseClass>
void MGDoFQuadAccessor<dim,BaseClass>::set_mg_dof_index (const unsigned int i,
							 const int index) const {
  Assert (dof_handler != 0, ExcInvalidObject());
  Assert (mg_dof_handler != 0, ExcInvalidObject());
				   // make sure a FE has been selected
				   // and enough room was reserved
  Assert (&dof_handler->get_fe() != 0, ExcInvalidObject());
  Assert (i<dof_handler->get_fe().dofs_per_quad,
	  ExcIndexRange (i, 0, dof_handler->get_fe().dofs_per_quad));

  mg_dof_handler->mg_levels[present_level]
    ->quad_dofs[present_index*dof_handler->get_fe().dofs_per_quad+i] = index;
};



template <int dim, typename BaseClass>
inline
int MGDoFQuadAccessor<dim,BaseClass>::mg_vertex_dof_index (const unsigned int vertex,
							   const unsigned int i) const {
  Assert (dof_handler != 0, ExcInvalidObject());
  Assert (mg_dof_handler != 0, ExcInvalidObject());
  Assert (&dof_handler->get_fe() != 0, ExcInvalidObject());
  Assert (vertex<4, ExcIndexRange (i,0,4));
  Assert (i<dof_handler->get_fe().dofs_per_vertex,
	  ExcIndexRange (i, 0, dof_handler->get_fe().dofs_per_vertex));
  
  return (mg_dof_handler->mg_vertex_dofs[vertex_index(vertex)]
	  .get_index (present_level, i, dof_handler->get_fe().dofs_per_vertex));
};


  
template <int dim, typename BaseClass>
void MGDoFQuadAccessor<dim,BaseClass>::set_mg_vertex_dof_index (const unsigned int vertex,
								const unsigned int i,
								const int index) const {
  Assert (dof_handler != 0, ExcInvalidObject());
  Assert (mg_dof_handler != 0, ExcInvalidObject());
  Assert (&dof_handler->get_fe() != 0, ExcInvalidObject());
  Assert (vertex<4, ExcIndexRange (i,0,4));
  Assert (i<dof_handler->get_fe().dofs_per_vertex,
	  ExcIndexRange (i, 0, dof_handler->get_fe().dofs_per_vertex));

  mg_dof_handler->mg_vertex_dofs[vertex_index(vertex)]
    .set_index (present_level, i, dof_handler->get_fe().dofs_per_vertex, index);
};



template <int dim, typename BaseClass>
void
MGDoFQuadAccessor<dim,BaseClass>::get_mg_dof_indices (vector<int> &dof_indices) const {
  Assert (dof_handler != 0, ExcInvalidObject());
  Assert (mg_dof_handler != 0, ExcInvalidObject());
  Assert (&dof_handler->get_fe() != 0, ExcInvalidObject());
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
      *next++ = mg_vertex_dof_index(vertex,d);
  for (unsigned int line=0; line<4; ++line)
    for (unsigned int d=0; d<dofs_per_line; ++d)
      *next++ = this->line(line)->mg_dof_index(d);
  for (unsigned int d=0; d<dofs_per_quad; ++d)
    *next++ = mg_dof_index(d);
};



template <int dim, typename BaseClass>
TriaIterator<dim,MGDoFLineAccessor<dim,LineAccessor<dim> > >
MGDoFQuadAccessor<dim,BaseClass>::line (const unsigned int i) const {
  Assert (i<4, ExcIndexRange (i, 0, 4));

  return TriaIterator<dim,MGDoFLineAccessor<dim,LineAccessor<dim> > >
    (
      tria,
      present_level,
      line_index (i),
      mg_dof_handler
    );
};



template <int dim, typename BaseClass>
TriaIterator<dim,MGDoFQuadAccessor<dim,BaseClass> >
MGDoFQuadAccessor<dim,BaseClass>::child (const unsigned int i) const {
  TriaIterator<dim,MGDoFQuadAccessor<dim,BaseClass> > q (tria,
							 present_level+1,
							 child_index (i),
							 mg_dof_handler);
  
#ifdef DEBUG
  if (q.state() != past_the_end)
    Assert (q->used(), typename TriaAccessor<dim>::ExcUnusedCellAsChild());
#endif
  return q;
};



template <int dim, typename BaseClass>
void
MGDoFQuadAccessor<dim,BaseClass>::copy_from (const MGDoFQuadAccessor<dim,BaseClass> &a) {
  DoFQuadAccessor<dim,BaseClass>::copy_from (a);
  set_mg_dof_handler (a.mg_dof_handler);
};




/* ------------------------ MGDoFHexAccessor --------------------------- */

template <int dim, typename BaseClass>
MGDoFHexAccessor<dim,BaseClass>::MGDoFHexAccessor (Triangulation<dim> *tria,
						   const int           level,
						   const int           index,
						   const AccessorData *local_data) :
		MGDoFAccessor<dim> (local_data),
		DoFHexAccessor<dim,BaseClass>(tria,level,index,local_data) {};



template <int dim, typename BaseClass>
inline
int MGDoFHexAccessor<dim,BaseClass>::mg_dof_index (const unsigned int i) const {
  Assert (dof_handler != 0, ExcInvalidObject());
  Assert (mg_dof_handler != 0, ExcInvalidObject());
				   // make sure a FE has been selected
				   // and enough room was reserved
  Assert (&dof_handler->get_fe() != 0, ExcInvalidObject());
  Assert (i<dof_handler->get_fe().dofs_per_hex,
	  ExcIndexRange (i, 0, dof_handler->get_fe().dofs_per_hex));

  return mg_dof_handler->mg_levels[present_level]
    ->hex_dofs[present_index*dof_handler->get_fe().dofs_per_hex+i];
};



template <int dim, typename BaseClass>
void MGDoFHexAccessor<dim,BaseClass>::set_mg_dof_index (const unsigned int i,
							const int index) const {
  Assert (dof_handler != 0, ExcInvalidObject());
  Assert (mg_dof_handler != 0, ExcInvalidObject());
				   // make sure a FE has been selected
				   // and enough room was reserved
  Assert (&dof_handler->get_fe() != 0, ExcInvalidObject());
  Assert (i<dof_handler->get_fe().dofs_per_hex,
	  ExcIndexRange (i, 0, dof_handler->get_fe().dofs_per_hex));

  mg_dof_handler->mg_levels[present_level]
    ->hex_dofs[present_index*dof_handler->get_fe().dofs_per_hex+i] = index;
};



template <int dim, typename BaseClass>
inline
int MGDoFHexAccessor<dim,BaseClass>::mg_vertex_dof_index (const unsigned int vertex,
							  const unsigned int i) const {
  Assert (dof_handler != 0, ExcInvalidObject());
  Assert (mg_dof_handler != 0, ExcInvalidObject());
  Assert (&dof_handler->get_fe() != 0, ExcInvalidObject());
  Assert (vertex<8, ExcIndexRange (i,0,8));
  Assert (i<dof_handler->get_fe().dofs_per_vertex,
	  ExcIndexRange (i, 0, dof_handler->get_fe().dofs_per_vertex));
  
  return (mg_dof_handler->mg_vertex_dofs[vertex_index(vertex)]
	  .get_index (present_level, i, dof_handler->get_fe().dofs_per_vertex));
};


  
template <int dim, typename BaseClass>
void MGDoFHexAccessor<dim,BaseClass>::set_mg_vertex_dof_index (const unsigned int vertex,
							       const unsigned int i,
							       const int index) const {
  Assert (dof_handler != 0, ExcInvalidObject());
  Assert (mg_dof_handler != 0, ExcInvalidObject());
  Assert (&dof_handler->get_fe() != 0, ExcInvalidObject());
  Assert (vertex<4, ExcIndexRange (i,0,4));
  Assert (i<dof_handler->get_fe().dofs_per_vertex,
	  ExcIndexRange (i, 0, dof_handler->get_fe().dofs_per_vertex));

  mg_dof_handler->mg_vertex_dofs[vertex_index(vertex)]
    .set_index (present_level, i, dof_handler->get_fe().dofs_per_vertex, index);
};



template <int dim, typename BaseClass>
void
MGDoFHexAccessor<dim,BaseClass>::get_mg_dof_indices (vector<int> &dof_indices) const {
  Assert (dof_handler != 0, ExcInvalidObject());
  Assert (mg_dof_handler != 0, ExcInvalidObject());
  Assert (&dof_handler->get_fe() != 0, ExcInvalidObject());
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
      *next++ = mg_vertex_dof_index(vertex,d);
  for (unsigned int line=0; line<12; ++line)
    for (unsigned int d=0; d<dofs_per_line; ++d)
      *next++ = this->line(line)->mg_dof_index(d);
  for (unsigned int quad=0; quad<12; ++quad)
    for (unsigned int d=0; d<dofs_per_quad; ++d)
      *next++ = this->quad(quad)->mg_dof_index(d);
  for (unsigned int d=0; d<dofs_per_hex; ++d)
    *next++ = mg_dof_index(d);
};



template <int dim, typename BaseClass>
TriaIterator<dim,MGDoFLineAccessor<dim,LineAccessor<dim> > >
MGDoFHexAccessor<dim,BaseClass>::line (const unsigned int i) const {
  Assert (i<12, ExcIndexRange (i, 0, 12));

  return TriaIterator<dim,MGDoFLineAccessor<dim,LineAccessor<dim> > >
    (
      tria,
      present_level,
      line_index (i),
      mg_dof_handler
    );
};



template <int dim, typename BaseClass>
TriaIterator<dim,MGDoFQuadAccessor<dim,QuadAccessor<dim> > >
MGDoFHexAccessor<dim,BaseClass>::quad (const unsigned int i) const {
  Assert (i<12, ExcIndexRange (i, 0, 6));

  return TriaIterator<dim,MGDoFQuadAccessor<dim,QuadAccessor<dim> > >
    (
      tria,
      present_level,
      quad_index (i),
      mg_dof_handler
    );
};



template <int dim, typename BaseClass>
TriaIterator<dim,MGDoFHexAccessor<dim,BaseClass> >
MGDoFHexAccessor<dim,BaseClass>::child (const unsigned int i) const {
  TriaIterator<dim,MGDoFHexAccessor<dim,BaseClass> > q (tria,
							present_level+1,
							child_index (i),
							mg_dof_handler);
  
#ifdef DEBUG
  if (q.state() != past_the_end)
    Assert (q->used(), typename TriaAccessor<dim>::ExcUnusedCellAsChild());
#endif
  return q;
};



template <int dim, typename BaseClass>
void
MGDoFHexAccessor<dim,BaseClass>::copy_from (const MGDoFHexAccessor<dim,BaseClass> &a) {
  DoFHexAccessor<dim,BaseClass>::copy_from (a);
  set_mg_dof_handler (a.mg_dof_handler);
};





/*------------------------- Functions: MGDoFCellAccessor -----------------------*/

template <int dim>
TriaIterator<dim,MGDoFCellAccessor<dim> >
MGDoFCellAccessor<dim>::neighbor (const unsigned int i) const {
  TriaIterator<dim,MGDoFCellAccessor<dim> > q (tria,
					       neighbor_level (i),
					       neighbor_index (i),
					       mg_dof_handler);
  
#ifdef DEBUG
  if (q.state() != past_the_end)
    Assert (q->used(), typename TriaAccessor<dim>::ExcUnusedCellAsNeighbor());
#endif
  return q;
};



template <int dim>
TriaIterator<dim,MGDoFCellAccessor<dim> >
MGDoFCellAccessor<dim>::child (const unsigned int i) const {
  TriaIterator<dim,MGDoFCellAccessor<dim> > q (tria,
					       present_level+1,
					       child_index (i),
					       mg_dof_handler);
  
#ifdef DEBUG
  if (q.state() != past_the_end)
    Assert (q->used(), typename TriaAccessor<dim>::ExcUnusedCellAsChild());
#endif
  return q;
};



#if deal_II_dimension == 1

template <>
MGDoFSubstructAccessor<1>::face_iterator
MGDoFCellAccessor<1>::face (const unsigned int) const {
  Assert (false, ExcNotUsefulForThisDimension());
  return 0;
};



template <>
void
MGDoFCellAccessor<1>::get_mg_dof_values (const Vector<double> &values,
					 Vector<double>       &dof_values) const {
  Assert (dof_handler != 0, ExcInvalidObject());
  Assert (mg_dof_handler != 0, ExcInvalidObject());
  Assert (&dof_handler->get_fe() != 0, ExcInvalidObject());
  Assert (dof_values.size() == dof_handler->get_fe().total_dofs,
	  ExcVectorDoesNotMatch());
  Assert (values.size() == dof_handler->n_dofs(),
	  ExcVectorDoesNotMatch());

  const unsigned int dofs_per_vertex = dof_handler->get_fe().dofs_per_vertex,
		     dofs_per_line   = dof_handler->get_fe().dofs_per_line;
  vector<double>::iterator next_dof_value=dof_values.begin();
  for (unsigned int vertex=0; vertex<2; ++vertex)
    for (unsigned int d=0; d<dofs_per_vertex; ++d)
      *next_dof_value++ = values(mg_vertex_dof_index(vertex,d));
  for (unsigned int d=0; d<dofs_per_line; ++d)
    *next_dof_value++ = values(mg_dof_index(d));
};

#endif



#if deal_II_dimension == 2

template <>
MGDoFSubstructAccessor<2>::face_iterator
MGDoFCellAccessor<2>::face (const unsigned int i) const {
  return line(i);
};



template <>
void
MGDoFCellAccessor<2>::get_mg_dof_values (const Vector<double> &values,
					 Vector<double>       &dof_values) const {
  Assert (dof_handler != 0, ExcInvalidObject());
  Assert (mg_dof_handler != 0, ExcInvalidObject());
  Assert (&dof_handler->get_fe() != 0, ExcInvalidObject());
  Assert (dof_values.size() == dof_handler->get_fe().total_dofs,
	  ExcVectorDoesNotMatch());
  Assert (values.size() == dof_handler->n_dofs(),
	  ExcVectorDoesNotMatch());

  const unsigned int dofs_per_vertex = dof_handler->get_fe().dofs_per_vertex,
		     dofs_per_line   = dof_handler->get_fe().dofs_per_line,
		     dofs_per_quad   = dof_handler->get_fe().dofs_per_quad;
  vector<double>::iterator next_dof_value=dof_values.begin();
  for (unsigned int vertex=0; vertex<4; ++vertex)
    for (unsigned int d=0; d<dofs_per_vertex; ++d)
      *next_dof_value++ = values(mg_vertex_dof_index(vertex,d));
  for (unsigned int line=0; line<4; ++line)
    for (unsigned int d=0; d<dofs_per_line; ++d)
      *next_dof_value++ = values(this->line(line)->mg_dof_index(d));
  for (unsigned int d=0; d<dofs_per_quad; ++d)
    *next_dof_value++ = values(mg_dof_index(d));
};

#endif





// explicit instantiations
#if deal_II_dimension == 1
template class MGDoFLineAccessor<1,CellAccessor<1> >;
template class MGDoFCellAccessor<1>;

template class TriaRawIterator<1,MGDoFCellAccessor<1> >;
template class TriaIterator<1,MGDoFCellAccessor<1> >;
template class TriaActiveIterator<1,MGDoFCellAccessor<1> >;
#endif

#if deal_II_dimension == 2
template class MGDoFLineAccessor<2,LineAccessor<2> >;
template class MGDoFQuadAccessor<2,QuadAccessor<2> >;
template class MGDoFQuadAccessor<2,CellAccessor<2> >;
template class MGDoFCellAccessor<2>;

template class TriaRawIterator<2,MGDoFLineAccessor<2,LineAccessor<2> > >;
template class TriaRawIterator<2,MGDoFCellAccessor<2> >;
template class TriaIterator<2,MGDoFLineAccessor<2,LineAccessor<2> > >;
template class TriaIterator<2,MGDoFCellAccessor<2> >;
template class TriaActiveIterator<2,MGDoFLineAccessor<2,LineAccessor<2> > >;
template class TriaActiveIterator<2,MGDoFCellAccessor<2> >;
#endif


#if deal_II_dimension == 3
template class MGDoFLineAccessor<3,LineAccessor<3> >;
template class MGDoFQuadAccessor<3,QuadAccessor<3> >;
template class MGDoFHexAccessor<3,HexAccessor<3> >;
template class MGDoFHexAccessor<3,CellAccessor<3> >;
template class MGDoFCellAccessor<3>;

template class TriaRawIterator<3,MGDoFLineAccessor<3,LineAccessor<3> > >;
template class TriaRawIterator<3,MGDoFQuadAccessor<3,QuadAccessor<3> > >;
template class TriaRawIterator<3,MGDoFHexAccessor<3,HexAccessor<3> > >;
template class TriaRawIterator<3,MGDoFCellAccessor<3> >;
template class TriaIterator<3,MGDoFLineAccessor<3,LineAccessor<3> > >;
template class TriaIterator<3,MGDoFQuadAccessor<3,QuadAccessor<3> > >;
template class TriaIterator<3,MGDoFCellAccessor<3> >;
template class TriaActiveIterator<3,MGDoFLineAccessor<3,LineAccessor<3> > >;
template class TriaActiveIterator<3,MGDoFQuadAccessor<3,QuadAccessor<3> > >;
template class TriaActiveIterator<3,MGDoFCellAccessor<3> >;
#endif


