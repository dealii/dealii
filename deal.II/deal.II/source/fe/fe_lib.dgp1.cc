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
template <> void FEDG_P1<deal_II_dimension>::initialize_matrices ();


#if deal_II_dimension == 1

template <>
FEDG_P1<1>::FEDG_P1 () :
		FEQ1Mapping<1> (0, 2, 0, 0, 1,
				vector<bool> (1, true))
{
//  initialize_matrices ();
};


template <>
void FEDG_P1<1>::initialize_matrices ()
{
				   // for restriction and prolongation matrices:
				   // note that we do not add up all the
				   // contributions since then we would get
				   // two summands per vertex in 1d (four
				   // in 2d, etc), but only one per line dof.
				   // We could accomplish for that by dividing
				   // the vertex dof values by 2 (4, etc), but
				   // would get into trouble at the boundary
				   // of the domain since there only one
				   // cell contributes to a vertex. Rather,
				   // we do not add up the contributions but
				   // set them right into the matrices!
  restriction[0](0,0) = 1.0;
  restriction[1](1,1) = 1.0;

  prolongation[0](0,0) = 1.0;
  prolongation[0](1,0) = 1./2.;
  prolongation[0](1,1) = 1./2.;

  prolongation[1](0,0) = 1./2.;
  prolongation[1](0,1) = 1./2.;
  prolongation[1](1,1) = 1.0;
};


template <>
double
FEDG_P1<1>::shape_value(const unsigned int i,
			 const Point<1>     &p) const
{
  Assert((i<dofs_per_cell), ExcIndexRange(i, 0, dofs_per_cell));
  switch (i)
    {
    case 0: return 1.-p(0);
    case 1: return p(0);
    }
  return 0.;
}


template <>
inline
Tensor<1,1>
FEDG_P1<1>::shape_grad(const unsigned int i,
			const Point<1>&) const
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
    }
  return Point<1>();
};


template <>
inline
Tensor<2,1>
FEDG_P1<1>::shape_grad_grad (const unsigned int i,
			     const Point<1> &) const
{
  Assert((i<dofs_per_cell), ExcIndexRange(i, 0, dofs_per_cell));
				   // second derivatives on the unit cell
				   // are always zero
  return Tensor<2,1>();
};


template <>
void FEDG_P1<1>::get_unit_support_points (vector<Point<1> >  &support_points) const
{
  FiniteElement<1>::get_unit_support_points (support_points);
};


template <>
void FEDG_P1<1>::get_support_points (const DoFHandler<1>::cell_iterator &cell,
				     vector<Point<1> >  &support_points) const
{
  FiniteElement<1>::get_support_points (cell, support_points);
};


template <>
void FEDG_P1<1>::get_face_support_points (const DoFHandler<1>::face_iterator &,
					  vector<Point<1> >  &) const
{
  Assert (false, ExcInternalError());
};


template <>
void FEDG_P1<1>::get_local_mass_matrix (const DoFHandler<1>::cell_iterator &cell,
					 FullMatrix<double> &local_mass_matrix) const
{
  Assert (local_mass_matrix.n() == dofs_per_cell,
	  ExcWrongFieldDimension(local_mass_matrix.n(),dofs_per_cell));
  Assert (local_mass_matrix.m() == dofs_per_cell,
	  ExcWrongFieldDimension(local_mass_matrix.m(),dofs_per_cell));

  const double h = cell->vertex(1)(0) - cell->vertex(0)(0);
  Assert (h>0, ExcJacobiDeterminantHasWrongSign());

  local_mass_matrix(0,0) = local_mass_matrix(1,1) = 1./3.*h;
  local_mass_matrix(0,1) = local_mass_matrix(1,0) = 1./6.*h;
};

#endif


#if deal_II_dimension == 2

template <>
FEDG_P1<2>::FEDG_P1 () :
		FEQ1Mapping<2> (0, 0, 3, 0, 1,
				vector<bool> (1, true))
{
//  initialize_matrices ();
};


template <>
void FEDG_P1<2>::initialize_matrices ()
{
  Assert(false, ExcNotImplemented());
};


template <>
inline
double
FEDG_P1<2>::shape_value (const unsigned int i,
			  const Point<2>& p) const
{
  Assert((i<dofs_per_cell), ExcIndexRange(i, 0, dofs_per_cell));
  switch (i)
    {
      case 0: return 1;
      case 1: return p(0);
      case 2: return p(1);
    }
  return 0.;
};


template <>
inline
Tensor<1,2>
FEDG_P1<2>::shape_grad (const unsigned int i,
			 const Point<2>&) const
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
    }
  return Point<2> ();
};


template <>
inline
Tensor<2,2>
FEDG_P1<2>::shape_grad_grad (const unsigned int i,
			      const Point<2> &) const
{
  Assert((i<dofs_per_cell), ExcIndexRange(i, 0, dofs_per_cell));
  return Tensor<2,2>();
};


template <>
void FEDG_P1<2>::get_local_mass_matrix (const DoFHandler<2>::cell_iterator &,
					 FullMatrix<double> &local_mass_matrix) const
{
  Assert(false, ExcNotImplemented ());
  Assert (local_mass_matrix.n() == dofs_per_cell,
	  ExcWrongFieldDimension(local_mass_matrix.n(),dofs_per_cell));
  Assert (local_mass_matrix.m() == dofs_per_cell,
	  ExcWrongFieldDimension(local_mass_matrix.m(),dofs_per_cell));
};


template <>
void FEDG_P1<2>::get_unit_support_points (vector<Point<2> > &unit_points) const
{
  Assert (unit_points.size() == dofs_per_cell,
	  ExcWrongFieldDimension (unit_points.size(), dofs_per_cell));

  unit_points[0] = Point<2> (.5,.5);
  unit_points[1] = Point<2> (1,0);
  unit_points[2] = Point<2> (0,1);
};


#endif


#if deal_II_dimension == 3

template <>
FEDG_P1<3>::FEDG_P1 () :
		FEQ1Mapping<3> (0, 0, 0, 4, 1,
				vector<bool> (1, true))
{
//  initialize_matrices ();
};


template <>
void FEDG_P1<3>::initialize_matrices ()
{
  Assert(false, ExcNotImplemented());
};


template <>
inline
double
FEDG_P1<3>::shape_value (const unsigned int i,
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
FEDG_P1<3>::shape_grad (const unsigned int i,
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
FEDG_P1<3>::shape_grad_grad (const unsigned int i,
			      const Point<3> &) const
{
  Assert((i<dofs_per_cell), ExcIndexRange(i, 0, dofs_per_cell));

  Tensor<2,3> return_value;
  return return_value;
};


template <>
void FEDG_P1<3>::get_local_mass_matrix (const DoFHandler<3>::cell_iterator &,
					 FullMatrix<double> &local_mass_matrix) const
{
  Assert (local_mass_matrix.n() == dofs_per_cell,
	  ExcWrongFieldDimension(local_mass_matrix.n(),dofs_per_cell));
  Assert (local_mass_matrix.m() == dofs_per_cell,
	  ExcWrongFieldDimension(local_mass_matrix.m(),dofs_per_cell));

  AssertThrow (false, ExcComputationNotUseful(3));
};


template <>
void FEDG_P1<3>::get_unit_support_points (vector<Point<3> > &unit_points) const
{
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
FEDG_P1<dim>::get_support_points (const typename DoFHandler<dim>::cell_iterator &cell,
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
FEDG_P1<dim>::get_face_support_points (const typename DoFHandler<dim>::face_iterator &face,
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

template class FEDG_P1<deal_II_dimension>;
