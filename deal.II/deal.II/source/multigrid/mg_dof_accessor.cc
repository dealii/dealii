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

template <int dim>
MGDoFObjectAccessor<1, dim>::MGDoFObjectAccessor (Triangulation<dim> *tria,
						     const int           level,
						     const int           index,
						     const AccessorData *local_data) :
		MGDoFAccessor<dim> (local_data),
		MGDoFObjectAccessor_Inheritance<1,dim>::BaseClass(tria,level,index,local_data) {};



template <int dim>
int MGDoFObjectAccessor<1, dim>::mg_dof_index (const unsigned int i) const {
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





template <int dim>
void MGDoFObjectAccessor<1, dim>::set_mg_dof_index (const unsigned int i,
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




template <int dim>
int MGDoFObjectAccessor<1, dim>::mg_vertex_dof_index (const unsigned int vertex,
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


  
template <int dim>
void MGDoFObjectAccessor<1, dim>::set_mg_vertex_dof_index (const unsigned int vertex,
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



template <int dim>
void
MGDoFObjectAccessor<1, dim>::get_mg_dof_indices (vector<int> &dof_indices) const {
  Assert (dof_handler != 0, ExcInvalidObject());
  Assert (mg_dof_handler != 0, ExcInvalidObject());
  Assert (&dof_handler->get_fe() != 0, ExcInvalidObject());
  Assert (dof_indices.size() == (2*dof_handler->get_fe().dofs_per_vertex +
				 dof_handler->get_fe().dofs_per_line),
	  DoFAccessor<dim>::ExcVectorDoesNotMatch());

  const unsigned int dofs_per_vertex = dof_handler->get_fe().dofs_per_vertex,
		     dofs_per_line   = dof_handler->get_fe().dofs_per_line;
  vector<int>::iterator next = dof_indices.begin();
  for (unsigned int vertex=0; vertex<2; ++vertex)
    for (unsigned int d=0; d<dofs_per_vertex; ++d)
      *next++ = mg_vertex_dof_index(vertex,d);
  for (unsigned int d=0; d<dofs_per_line; ++d)
    *next++ = mg_dof_index(d);
};



template <int dim>
TriaIterator<dim,MGDoFObjectAccessor<1, dim> >
MGDoFObjectAccessor<1, dim>::child (const unsigned int i) const {
  TriaIterator<dim,MGDoFObjectAccessor<1, dim> > q (tria,
							 present_level+1,
							 child_index (i),
							 mg_dof_handler);
  
#ifdef DEBUG
  if (q.state() != past_the_end)
    Assert (q->used(), typename TriaAccessor<dim>::ExcUnusedCellAsChild());
#endif
  return q;
};



template <int dim>
void
MGDoFObjectAccessor<1, dim>::copy_from (const MGDoFObjectAccessor<1, dim> &a)
{
  DoFObjectAccessor<1, dim>::copy_from (a);
  set_mg_dof_handler (a.mg_dof_handler);
};




/* ------------------------ MGDoFQuadAccessor --------------------------- */

template <int dim>
MGDoFObjectAccessor<2, dim>::MGDoFObjectAccessor (Triangulation<dim> *tria,
						  const int           level,
						  const int           index,
						  const AccessorData *local_data) :
		MGDoFAccessor<dim> (local_data),
		MGDoFObjectAccessor_Inheritance<2,dim>::BaseClass(tria,level,index,local_data)
{};



template <int dim>
inline
int MGDoFObjectAccessor<2, dim>::mg_dof_index (const unsigned int i) const {
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



template <int dim>
void MGDoFObjectAccessor<2, dim>::set_mg_dof_index (const unsigned int i,
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



template <int dim>
inline
int MGDoFObjectAccessor<2, dim>::mg_vertex_dof_index (const unsigned int vertex,
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


  
template <int dim>
void MGDoFObjectAccessor<2, dim>::set_mg_vertex_dof_index (const unsigned int vertex,
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



template <int dim>
void
MGDoFObjectAccessor<2, dim>::get_mg_dof_indices (vector<int> &dof_indices) const {
  Assert (dof_handler != 0, ExcInvalidObject());
  Assert (mg_dof_handler != 0, ExcInvalidObject());
  Assert (&dof_handler->get_fe() != 0, ExcInvalidObject());
  Assert (dof_indices.size() == (4*dof_handler->get_fe().dofs_per_vertex +
				 4*dof_handler->get_fe().dofs_per_line +
				 dof_handler->get_fe().dofs_per_quad),
	  DoFAccessor<2>::ExcVectorDoesNotMatch());

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



template <int dim>
TriaIterator<dim,MGDoFObjectAccessor<1, dim> >
MGDoFObjectAccessor<2, dim>::line (const unsigned int i) const {
  Assert (i<4, ExcIndexRange (i, 0, 4));

  return TriaIterator<dim,MGDoFObjectAccessor<1, dim> >
    (
      tria,
      present_level,
      line_index (i),
      mg_dof_handler
    );
};



template <int dim>
TriaIterator<dim,MGDoFObjectAccessor<2, dim> >
MGDoFObjectAccessor<2, dim>::child (const unsigned int i) const {
  TriaIterator<dim,MGDoFObjectAccessor<2, dim> > q (tria,
							 present_level+1,
							 child_index (i),
							 mg_dof_handler);
  
#ifdef DEBUG
  if (q.state() != past_the_end)
    Assert (q->used(), typename TriaAccessor<dim>::ExcUnusedCellAsChild());
#endif
  return q;
};



template <int dim>
void
MGDoFObjectAccessor<2, dim>::copy_from (const MGDoFObjectAccessor<2, dim> &a) {
  DoFObjectAccessor<2, dim>::copy_from (a);
  set_mg_dof_handler (a.mg_dof_handler);
};




/* ------------------------ MGDoFHexAccessor --------------------------- */

template <int dim>
MGDoFObjectAccessor<3, dim>::MGDoFObjectAccessor (Triangulation<dim> *tria,
						   const int           level,
						   const int           index,
						   const AccessorData *local_data) :
		MGDoFAccessor<dim> (local_data),
		MGDoFObjectAccessor_Inheritance<3,dim>::BaseClass(tria,level,index,local_data) {};



template <int dim>
inline
int MGDoFObjectAccessor<3, dim>::mg_dof_index (const unsigned int i) const {
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



template <int dim>
void MGDoFObjectAccessor<3, dim>::set_mg_dof_index (const unsigned int i,
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



template <int dim>
inline
int MGDoFObjectAccessor<3, dim>::mg_vertex_dof_index (const unsigned int vertex,
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


  
template <int dim>
void MGDoFObjectAccessor<3, dim>::set_mg_vertex_dof_index (const unsigned int vertex,
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



template <int dim>
void
MGDoFObjectAccessor<3, dim>::get_mg_dof_indices (vector<int> &dof_indices) const {
  Assert (dof_handler != 0, ExcInvalidObject());
  Assert (mg_dof_handler != 0, ExcInvalidObject());
  Assert (&dof_handler->get_fe() != 0, ExcInvalidObject());
  Assert (dof_indices.size() == (8*dof_handler->get_fe().dofs_per_vertex +
				 12*dof_handler->get_fe().dofs_per_line +
				 6*dof_handler->get_fe().dofs_per_quad +
				 dof_handler->get_fe().dofs_per_hex),
	  DoFAccessor<3>::ExcVectorDoesNotMatch());

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



template <int dim>
TriaIterator<dim,MGDoFObjectAccessor<1, dim> >
MGDoFObjectAccessor<3, dim>::line (const unsigned int i) const {
  Assert (i<12, ExcIndexRange (i, 0, 12));

  return TriaIterator<dim,MGDoFObjectAccessor<1, dim> >
    (
      tria,
      present_level,
      line_index (i),
      mg_dof_handler
    );
};



template <int dim>
TriaIterator<dim,MGDoFObjectAccessor<2, dim> >
MGDoFObjectAccessor<3, dim>::quad (const unsigned int i) const {
  Assert (i<12, ExcIndexRange (i, 0, 6));

  return TriaIterator<dim,MGDoFObjectAccessor<2, dim> >
    (
      tria,
      present_level,
      quad_index (i),
      mg_dof_handler
    );
};



template <int dim>
TriaIterator<dim,MGDoFObjectAccessor<3, dim> >
MGDoFObjectAccessor<3, dim>::child (const unsigned int i) const {
  TriaIterator<dim,MGDoFObjectAccessor<3, dim> > q (tria,
						    present_level+1,
						    child_index (i),
						    mg_dof_handler);
  
#ifdef DEBUG
  if (q.state() != past_the_end)
    Assert (q->used(), typename TriaAccessor<dim>::ExcUnusedCellAsChild());
#endif
  return q;
};



template <int dim>
void
MGDoFObjectAccessor<3, dim>::copy_from (const MGDoFObjectAccessor<3, dim> &a) {
  DoFObjectAccessor<3, dim>::copy_from (a);
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
MGDoFCellAccessor<1>::face_iterator
MGDoFCellAccessor<1>::face (const unsigned int) const
{
  Assert (false, ExcNotUsefulForThisDimension());
  return face_iterator();
};



template <>
void
MGDoFCellAccessor<1>::get_mg_dof_values (const Vector<double> &values,
					 Vector<double>       &dof_values) const {
  Assert (dof_handler != 0, DoFAccessor<1>::ExcInvalidObject());
  Assert (mg_dof_handler != 0, DoFAccessor<1>::ExcInvalidObject());
  Assert (&dof_handler->get_fe() != 0, DoFAccessor<1>::ExcInvalidObject());
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
MGDoFCellAccessor<2>::face_iterator
MGDoFCellAccessor<2>::face (const unsigned int i) const
{
  return line(i);
};



template <>
void
MGDoFCellAccessor<2>::get_mg_dof_values (const Vector<double> &values,
					 Vector<double>       &dof_values) const {
  Assert (dof_handler != 0, DoFAccessor<2>::ExcInvalidObject());
  Assert (mg_dof_handler != 0, DoFAccessor<2>::ExcInvalidObject());
  Assert (&dof_handler->get_fe() != 0, DoFAccessor<2>::ExcInvalidObject());
  Assert (dof_values.size() == dof_handler->get_fe().total_dofs,
	  DoFAccessor<2>::ExcVectorDoesNotMatch());
  Assert (values.size() == mg_dof_handler->n_dofs(present_level),
	  DoFAccessor<2>::ExcVectorDoesNotMatch());

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
template class MGDoFObjectAccessor<1, 1>;
template class MGDoFCellAccessor<1>;

template class TriaRawIterator<1,MGDoFCellAccessor<1> >;
template class TriaIterator<1,MGDoFCellAccessor<1> >;
template class TriaActiveIterator<1,MGDoFCellAccessor<1> >;
#endif

#if deal_II_dimension == 2
template class MGDoFObjectAccessor<1, 2>;
template class MGDoFObjectAccessor<2, 2>;
template class MGDoFCellAccessor<2>;

template class TriaRawIterator   <2,MGDoFObjectAccessor<1, 2> >;
template class TriaIterator      <2,MGDoFObjectAccessor<1, 2> >;
template class TriaActiveIterator<2,MGDoFObjectAccessor<1, 2> >;

template class TriaRawIterator   <2,MGDoFCellAccessor<2> >;
template class TriaIterator      <2,MGDoFCellAccessor<2> >;
template class TriaActiveIterator<2,MGDoFCellAccessor<2> >;
#endif


#if deal_II_dimension == 3
template class MGDoFObjectAccessor<1, 3>;
template class MGDoFObjectAccessor<2, 3>;
template class MGDoFObjectAccessor<3, 3>;
template class MGDoFCellAccessor<3>;

template class TriaRawIterator   <3,MGDoFObjectAccessor<1, 3> >;
template class TriaIterator      <3,MGDoFObjectAccessor<1, 3> >;
template class TriaActiveIterator<3,MGDoFObjectAccessor<1, 3> >;

template class TriaRawIterator   <3,MGDoFObjectAccessor<2, 3> >;
template class TriaIterator      <3,MGDoFObjectAccessor<2, 3> >;
template class TriaActiveIterator<3,MGDoFObjectAccessor<2, 3> >;

template class TriaRawIterator   <3,MGDoFCellAccessor<3> >;
template class TriaIterator      <3,MGDoFCellAccessor<3> >;
template class TriaActiveIterator<3,MGDoFCellAccessor<3> >;
#endif


