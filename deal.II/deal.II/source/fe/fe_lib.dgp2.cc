//----------------------------  $RCSFile$  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  $RCSFile$  ---------------------------


#include <fe/fe_lib.dgp.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_accessor.h>
#include <grid/geometry_info.h>
#include <algorithm>


// declare explicit specializations before use:
template <> void FEDG_P2<deal_II_dimension>::initialize_matrices ();


#if deal_II_dimension == 1

template <>
FEDG_P2<1>::FEDG_P2 () :
		FEQ1Mapping<1> (0, 3, 0, 0, 1,
				vector<bool> (1, true))
{
//  initialize_matrices ();
};


template <>
void FEDG_P2<1>::initialize_matrices ()
{
  Assert(false, ExcNotImplemented());
};


template <>
double
FEDG_P2<1>::shape_value(const unsigned int i,
			 const Point<1>     &p) const
{
  Assert((i<dofs_per_cell), ExcIndexRange(i, 0, dofs_per_cell));
  switch (i)
    {
      case 0: return 1.;
      case 1: return p(0);
      case 2: return p(0)*p(0);
    }
  return 0.;
}


template <>
inline
Tensor<1,1>
FEDG_P2<1>::shape_grad(const unsigned int i,
			const Point<1>&p) const
{
  Assert((i<dofs_per_cell), ExcIndexRange(i, 0, dofs_per_cell));
				   // originally, the return type of the
				   // function was Point<dim>, so we
				   // still construct it as that. it should
				   // make no difference in practice,
				   // however
  switch (i)
    {
      case 0: return Point<1>(-1.);
      case 1: return Point<1>(1.);
      case 2: return Point<1>(2.*p(0));
	    
    }
  return Point<1>();
};


template <>
inline
Tensor<2,1>
FEDG_P2<1>::shape_grad_grad (const unsigned int i,
			     const Point<1> &) const
{
  Assert(false, ExcNotImplemented());
  Assert((i<dofs_per_cell), ExcIndexRange(i, 0, dofs_per_cell));
  return Tensor<2,1>();
};


template <>
void FEDG_P2<1>::get_unit_support_points (vector<Point<1> >  &support_points) const
{
  FiniteElement<1>::get_unit_support_points (support_points);
};


template <>
void FEDG_P2<1>::get_support_points (const DoFHandler<1>::cell_iterator &cell,
				     vector<Point<1> >  &support_points) const
{
  FiniteElement<1>::get_support_points (cell, support_points);
};


template <>
void FEDG_P2<1>::get_face_support_points (const DoFHandler<1>::face_iterator &,
					  vector<Point<1> >  &) const
{
  Assert (false, ExcInternalError());
};


template <>
void FEDG_P2<1>::get_local_mass_matrix (const DoFHandler<1>::cell_iterator &,
					 FullMatrix<double> &local_mass_matrix) const
{
  Assert(false, ExcNotImplemented());
  Assert (local_mass_matrix.n() == dofs_per_cell,
	  ExcWrongFieldDimension(local_mass_matrix.n(),dofs_per_cell));
  Assert (local_mass_matrix.m() == dofs_per_cell,
	  ExcWrongFieldDimension(local_mass_matrix.m(),dofs_per_cell));
};

#endif


#if deal_II_dimension == 2

template <>
FEDG_P2<2>::FEDG_P2 () :
		FEQ1Mapping<2> (0, 0, 6, 0, 1,
				vector<bool> (1, true))
{
//  initialize_matrices ();
};


template <>
void FEDG_P2<2>::initialize_matrices ()
{
  Assert(false, ExcNotImplemented());
};


template <>
inline
double
FEDG_P2<2>::shape_value (const unsigned int i,
			  const Point<2>& p) const
{
  Assert((i<dofs_per_cell), ExcIndexRange(i, 0, dofs_per_cell));
  switch (i)
    {
      case 0: return 1;
      case 1: return p(0);
      case 2: return p(1);
      case 3: return p(0)*p(0);
      case 4: return p(0)*p(1);
      case 5: return p(1)*p(1);
    }
  return 0.;
};


template <>
inline
Tensor<1,2>
FEDG_P2<2>::shape_grad (const unsigned int i,
			 const Point<2>& p) const
{
  Assert((i<dofs_per_cell), ExcIndexRange(i, 0, dofs_per_cell));
				   // originally, the return type of the
				   // function was Point<dim>, so we
				   // still construct it as that. it should
				   // make no difference in practice,
				   // however
  switch (i)
    {
      case 0: return Point<2> (0,0);
      case 1: return Point<2> (1,0);
      case 2: return Point<2> (0,1);
      case 3: return Point<2> (2*p(0),0);
      case 4: return Point<2> (p(1),p(0));
      case 5: return Point<2> (0,2*p(1));
    }
  return Point<2> ();
};


template <>
inline
Tensor<2,2>
FEDG_P2<2>::shape_grad_grad (const unsigned int i,
			      const Point<2> &) const
{
  Assert(false, ExcNotImplemented());
  Assert((i<dofs_per_cell), ExcIndexRange(i, 0, dofs_per_cell));
  return Tensor<2,2>();
};


template <>
void FEDG_P2<2>::get_local_mass_matrix (const DoFHandler<2>::cell_iterator &,
					 FullMatrix<double> &) const
{
  Assert(false, ExcNotImplemented ());
};


template <>
void FEDG_P2<2>::get_unit_support_points (vector<Point<2> > &unit_points) const
{
  Assert (unit_points.size() == dofs_per_cell,
	  ExcWrongFieldDimension (unit_points.size(), dofs_per_cell));

  unit_points[0] = Point<2> (.5,.5);
  unit_points[1] = Point<2> (1,0);
  unit_points[2] = Point<2> (0,1);
  unit_points[3] = Point<2> (1,0);
  unit_points[4] = Point<2> (0,1);
  unit_points[5] = Point<2> (1,1);
};


#endif


#if deal_II_dimension == 3

template <>
FEDG_P2<3>::FEDG_P2 () :
		FEQ1Mapping<3> (0, 0, 0, 4, 1,
				vector<bool> (1, true))
{
  Assert(false, ExcNotImplemented ());
//  initialize_matrices ();
};


template <>
void FEDG_P2<3>::initialize_matrices ()
{
  Assert(false, ExcNotImplemented());
};


template <>
inline
double
FEDG_P2<3>::shape_value (const unsigned int i,
			 const Point<3>& p) const
{
  Assert((i<dofs_per_cell), ExcIndexRange(i, 0, dofs_per_cell));
  switch (i)
    {
      case 0: return 1.;
      case 1: return p(0);
      case 2: return p(1);
      case 3: return p(2);
    }
  return 0.;
};


template <>
inline
Tensor<1,3>
FEDG_P2<3>::shape_grad (const unsigned int i,
			 const Point<3>&) const
{
  Assert((i<dofs_per_cell), ExcIndexRange(i, 0, dofs_per_cell));
				   // originally, the return type of the
				   // function was Point<dim>, so we
				   // still construct it as that. it should
				   // make no difference in practice,
				   // however
  switch (i)
    {
      case 0: return Point<3>(0,0,0);
      case 1: return Point<3>(1,0,0);
      case 2: return Point<3>(0,1,0);
      case 3: return Point<3>(0,0,1);
    }
  return Point<3> ();
};


template <>
inline
Tensor<2,3>
FEDG_P2<3>::shape_grad_grad (const unsigned int i,
			      const Point<3> &) const
{
  Assert((i<dofs_per_cell), ExcIndexRange(i, 0, dofs_per_cell));

  Tensor<2,3> return_value;
  return return_value;
};


template <>
void FEDG_P2<3>::get_local_mass_matrix (const DoFHandler<3>::cell_iterator &,
					 FullMatrix<double> &local_mass_matrix) const
{
  Assert (local_mass_matrix.n() == dofs_per_cell,
	  ExcWrongFieldDimension(local_mass_matrix.n(),dofs_per_cell));
  Assert (local_mass_matrix.m() == dofs_per_cell,
	  ExcWrongFieldDimension(local_mass_matrix.m(),dofs_per_cell));

  throw ExcComputationNotUseful(3);
};


template <>
void FEDG_P2<3>::get_unit_support_points (vector<Point<3> > &unit_points) const {
  Assert (unit_points.size() == dofs_per_cell,
	  ExcWrongFieldDimension (unit_points.size(), dofs_per_cell));

  unit_points[0] = Point<3> (.5,.5,.5);
  unit_points[1] = Point<3> (1,0,0);
  unit_points[2] = Point<3> (0,1,0);
  unit_points[3] = Point<3> (0,0,1);
};


#endif


template <int dim>
void
FEDG_P2<dim>::get_support_points (const typename DoFHandler<dim>::cell_iterator &cell,
				   vector<Point<dim> >  &support_points) const
{
  Assert (support_points.size() == dofs_per_cell,
	  typename FiniteElementBase<dim>::ExcWrongFieldDimension (support_points.size(),
								   dofs_per_cell));
  
  for (unsigned int vertex=0; vertex<GeometryInfo<dim>::vertices_per_cell; ++vertex)
    support_points[vertex] = cell->vertex(vertex);
};


template <int dim>
void
FEDG_P2<dim>::get_face_support_points (const typename DoFHandler<dim>::face_iterator &face,
				       vector<Point<dim> >  &support_points) const
{
  Assert ((support_points.size() == dofs_per_face) &&
	  (support_points.size() == GeometryInfo<dim>::vertices_per_face),
	  typename FiniteElementBase<dim>::ExcWrongFieldDimension (support_points.size(),
								   GeometryInfo<dim>::vertices_per_face));

  for (unsigned int vertex=0; vertex<dofs_per_face; ++vertex)
    support_points[vertex] = face->vertex(vertex);
};


// explicit instantiations

template class FEDG_P2<deal_II_dimension>;
