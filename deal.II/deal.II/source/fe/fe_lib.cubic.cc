/* $Id$ */
/* Copyright W. Bangerth, University of Heidelberg, 1998 */


#include <fe/fe_lib.lagrange.h>
#include <grid/tria_iterator.h>
#include <grid/dof_accessor.h>
#include <grid/geometry_info.h>
#include <algorithm>






#if deal_II_dimension == 1

template <>
FECubicSub<1>::FECubicSub () :
		FELinearMapping<1> (1, 2) {
  prolongation[0](0,0) = 1.0;
  prolongation[0](0,1) = 0.0;
  prolongation[0](0,2) = 0.0;
  prolongation[0](0,3) = 0.0;
  prolongation[0](1,0) = -1.0/16.0;
  prolongation[0](1,1) = -1.0/16.0;
  prolongation[0](1,2) = 9.0/16.0;
  prolongation[0](1,3) = 9.0/16.0;
  prolongation[0](2,0) = 5.0/16.0;
  prolongation[0](2,1) = 1.0/16.0;
  prolongation[0](2,2) = 15.0/16.0;
  prolongation[0](2,3) = -5.0/16.0;
  prolongation[0](3,0) = 0.0;
  prolongation[0](3,1) = 0.0;
  prolongation[0](3,2) = 1.0;
  prolongation[0](3,3) = 0.0;
  prolongation[1](0,0) = -1.0/16.0;
  prolongation[1](0,1) = -1.0/16.0;
  prolongation[1](0,2) = 9.0/16.0;
  prolongation[1](0,3) = 9.0/16.0;
  prolongation[1](1,0) = 0.0;
  prolongation[1](1,1) = 1.0;
  prolongation[1](1,2) = 0.0;
  prolongation[1](1,3) = 0.0;
  prolongation[1](2,0) = 0.0;
  prolongation[1](2,1) = 0.0;
  prolongation[1](2,2) = 0.0;
  prolongation[1](2,3) = 1.0;
  prolongation[1](3,0) = 1.0/16.0;
  prolongation[1](3,1) = 5.0/16.0;
  prolongation[1](3,2) = -5.0/16.0;
  prolongation[1](3,3) = 15.0/16.0;

  restriction[0](0,0) = 1.0;
  restriction[0](0,1) = 0.0;
  restriction[0](0,2) = 0.0;
  restriction[0](0,3) = 0.0;
  restriction[0](1,0) = 0.0;
  restriction[0](1,1) = 0.0;
  restriction[0](1,2) = 0.0;
  restriction[0](1,3) = 0.0;
  restriction[0](2,0) = 0.0;
  restriction[0](2,1) = 0.0;
  restriction[0](2,2) = 0.0;
  restriction[0](2,3) = 1.0;
  restriction[0](3,0) = 0.0;
  restriction[0](3,1) = 0.0;
  restriction[0](3,2) = 0.0;
  restriction[0](3,3) = 0.0;
  restriction[1](0,0) = 0.0;
  restriction[1](0,1) = 0.0;
  restriction[1](0,2) = 0.0;
  restriction[1](0,3) = 0.0;
  restriction[1](1,0) = 0.0;
  restriction[1](1,1) = 1.0;
  restriction[1](1,2) = 0.0;
  restriction[1](1,3) = 0.0;
  restriction[1](2,0) = 0.0;
  restriction[1](2,1) = 0.0;
  restriction[1](2,2) = 0.0;
  restriction[1](2,3) = 0.0;
  restriction[1](3,0) = 0.0;
  restriction[1](3,1) = 0.0;
  restriction[1](3,2) = 1.0;
  restriction[1](3,3) = 0.0;
};



template <>
double
FECubicSub<1>::shape_value (const unsigned int i,
			    const Point<1>     &p) const
{
  Assert((i<total_dofs), ExcInvalidIndex(i));
  const double xi = p(0);
  switch (i)
    {
      case 0: return -9.0/2.0*xi*xi*xi+9.0*xi*xi-11.0/2.0*xi+1.0;
      case 1: return 9.0/2.0*xi*xi*xi-9.0/2.0*xi*xi+xi;
      case 2: return 27.0/2.0*xi*xi*xi-45.0/2.0*xi*xi+9.0*xi;
      case 3: return -27.0/2.0*xi*xi*xi+18.0*xi*xi-9.0/2.0*xi;
    }
  return 0.;
};



template <>
Point<1>
FECubicSub<1>::shape_grad (const unsigned int i,
			   const Point<1>    &p) const
{
  Assert((i<total_dofs), ExcInvalidIndex(i));
  const double xi = p(0);
  switch (i)
    {
      case 0: return Point<1>(-27.0/2.0*xi*xi+18.0*xi-11.0/2.0);
      case 1: return Point<1>(27.0/2.0*xi*xi-9.0*xi+1.0);
      case 2: return Point<1>(81.0/2.0*xi*xi-45.0*xi+9.0);
      case 3: return Point<1>(-81.0/2.0*xi*xi+36.0*xi-9.0/2.0);
    }
  return Point<1>();
};



template <>
Tensor<2,1>
FECubicSub<1>::shape_grad_grad (const unsigned int i,
				const Point<1>    &p) const
{
  Assert (i<total_dofs, ExcInvalidIndex(i));

  const double xi = p(0);
  Tensor<2,1> return_value;
  switch (i) 
    {
      case 0: return_value[0][0] = -27.0*xi+18.0;
      case 1: return_value[0][0] = 27.0*xi-9.0;
      case 2: return_value[0][0] = 81.0*xi-45.0;
      case 3: return_value[0][0] = -81.0*xi+36.0;
    };

  return return_value;
};



template <>
void FECubicSub<1>::get_unit_support_points (vector<Point<1> > &unit_points) const {
  FiniteElement<1>::get_unit_support_points (unit_points);
};



template <>
void FECubicSub<1>::get_support_points (const typename DoFHandler<1>::cell_iterator &cell,
					   const Boundary<1>  &boundary,
					   vector<Point<1> >  &support_points) const {
  FiniteElement<1>::get_support_points (cell, boundary, support_points);
};



template <>
void FECubicSub<1>::get_face_support_points (const typename DoFHandler<1>::face_iterator &,
					     const Boundary<1>  &,
					     vector<Point<1> >  &) const {
  Assert (false, ExcInternalError());
};



template <>
void FECubicSub<1>::get_local_mass_matrix (const DoFHandler<1>::cell_iterator &cell,
					       const Boundary<1> &,
					       dFMatrix &local_mass_matrix) const {
  Assert (local_mass_matrix.n() == total_dofs,
	  ExcWrongFieldDimension(local_mass_matrix.n(),total_dofs));
  Assert (local_mass_matrix.m() == total_dofs,
	  ExcWrongFieldDimension(local_mass_matrix.m(),total_dofs));

  const double h = cell->vertex(1)(0) - cell->vertex(0)(0);
  Assert (h>0, ExcJacobiDeterminantHasWrongSign());

  const double t1 = 8.0/105.0*h;
  const double t2 = 19.0/1680.0*h;
  const double t3 = 33.0/560.0*h;
  const double t4 = 3.0/140.0*h;
  const double t5 = 27.0/70.0*h;
  const double t6 = 27.0/560.0*h;
  local_mass_matrix(0,0) = t1;
  local_mass_matrix(0,1) = t2;
  local_mass_matrix(0,2) = t3;
  local_mass_matrix(0,3) = -t4;
  local_mass_matrix(1,0) = t2;
  local_mass_matrix(1,1) = t1;
  local_mass_matrix(1,2) = -t4;
  local_mass_matrix(1,3) = t3;
  local_mass_matrix(2,0) = t3;
  local_mass_matrix(2,1) = -t4;
  local_mass_matrix(2,2) = t5;
  local_mass_matrix(2,3) = -t6;
  local_mass_matrix(3,0) = -t4;
  local_mass_matrix(3,1) = t3;
  local_mass_matrix(3,2) = -t6;
  local_mass_matrix(3,3) = t5;
};

#endif




#if deal_II_dimension == 2

template <>
FECubicSub<2>::FECubicSub () :
		FELinearMapping<2> (1, 2, 4)
{
  interface_constraints(0,0) = -1.0/16.0;
  interface_constraints(0,1) = -1.0/16.0;
  interface_constraints(0,2) = 9.0/16.0;
  interface_constraints(0,3) = 9.0/16.0;
  interface_constraints(1,0) = 5.0/16.0;
  interface_constraints(1,1) = 1.0/16.0;
  interface_constraints(1,2) = 15.0/16.0;
  interface_constraints(1,3) = -5.0/16.0;
  interface_constraints(2,2) = 1.0;
  interface_constraints(3,3) = 1.0;
  interface_constraints(4,0) = 1.0/16.0;
  interface_constraints(4,1) = 5.0/16.0;
  interface_constraints(4,2) = -5.0/16.0;
  interface_constraints(4,3) = 15.0/16.0;

  prolongation[0](0,0) = 1.0;
  prolongation[0](1,0) = -1.0/16.0;
  prolongation[0](1,1) = -1.0/16.0;
  prolongation[0](1,4) = 9.0/16.0;
  prolongation[0](1,5) = 9.0/16.0;
  prolongation[0](2,0) = 1.0/256.0;
  prolongation[0](2,1) = 1.0/256.0;
  prolongation[0](2,2) = 1.0/256.0;
  prolongation[0](2,3) = 1.0/256.0;
  prolongation[0](2,4) = -9.0/256.0;
  prolongation[0](2,5) = -9.0/256.0;
  prolongation[0](2,6) = -9.0/256.0;
  prolongation[0](2,7) = -9.0/256.0;
  prolongation[0](2,8) = -9.0/256.0;
  prolongation[0](2,9) = -9.0/256.0;
  prolongation[0](2,10) = -9.0/256.0;
  prolongation[0](2,11) = -9.0/256.0;
  prolongation[0](2,12) = 81.0/256.0;
  prolongation[0](2,13) = 81.0/256.0;
  prolongation[0](2,14) = 81.0/256.0;
  prolongation[0](2,15) = 81.0/256.0;
  prolongation[0](3,0) = -1.0/16.0;
  prolongation[0](3,3) = -1.0/16.0;
  prolongation[0](3,10) = 9.0/16.0;
  prolongation[0](3,11) = 9.0/16.0;
  prolongation[0](4,0) = 5.0/16.0;
  prolongation[0](4,1) = 1.0/16.0;
  prolongation[0](4,4) = 15.0/16.0;
  prolongation[0](4,5) = -5.0/16.0;
  prolongation[0](5,4) = 1.0;
  prolongation[0](6,0) = -5.0/256.0;
  prolongation[0](6,1) = -5.0/256.0;
  prolongation[0](6,2) = -1.0/256.0;
  prolongation[0](6,3) = -1.0/256.0;
  prolongation[0](6,4) = 45.0/256.0;
  prolongation[0](6,5) = 45.0/256.0;
  prolongation[0](6,6) = -15.0/256.0;
  prolongation[0](6,7) = 5.0/256.0;
  prolongation[0](6,8) = 9.0/256.0;
  prolongation[0](6,9) = 9.0/256.0;
  prolongation[0](6,10) = -15.0/256.0;
  prolongation[0](6,11) = 5.0/256.0;
  prolongation[0](6,12) = 135.0/256.0;
  prolongation[0](6,13) = 135.0/256.0;
  prolongation[0](6,14) = -45.0/256.0;
  prolongation[0](6,15) = -45.0/256.0;
  prolongation[0](7,6) = -1.0/16.0;
  prolongation[0](7,10) = -1.0/16.0;
  prolongation[0](7,12) = 9.0/16.0;
  prolongation[0](7,13) = 9.0/16.0;
  prolongation[0](8,0) = -5.0/256.0;
  prolongation[0](8,1) = -1.0/256.0;
  prolongation[0](8,2) = -1.0/256.0;
  prolongation[0](8,3) = -5.0/256.0;
  prolongation[0](8,4) = -15.0/256.0;
  prolongation[0](8,5) = 5.0/256.0;
  prolongation[0](8,6) = 9.0/256.0;
  prolongation[0](8,7) = 9.0/256.0;
  prolongation[0](8,8) = -15.0/256.0;
  prolongation[0](8,9) = 5.0/256.0;
  prolongation[0](8,10) = 45.0/256.0;
  prolongation[0](8,11) = 45.0/256.0;
  prolongation[0](8,12) = 135.0/256.0;
  prolongation[0](8,13) = -45.0/256.0;
  prolongation[0](8,14) = -45.0/256.0;
  prolongation[0](8,15) = 135.0/256.0;
  prolongation[0](9,4) = -1.0/16.0;
  prolongation[0](9,8) = -1.0/16.0;
  prolongation[0](9,12) = 9.0/16.0;
  prolongation[0](9,15) = 9.0/16.0;
  prolongation[0](10,0) = 5.0/16.0;
  prolongation[0](10,3) = 1.0/16.0;
  prolongation[0](10,10) = 15.0/16.0;
  prolongation[0](10,11) = -5.0/16.0;
  prolongation[0](11,10) = 1.0;
  prolongation[0](12,0) = 25.0/256.0;
  prolongation[0](12,1) = 5.0/256.0;
  prolongation[0](12,2) = 1.0/256.0;
  prolongation[0](12,3) = 5.0/256.0;
  prolongation[0](12,4) = 75.0/256.0;
  prolongation[0](12,5) = -25.0/256.0;
  prolongation[0](12,6) = 15.0/256.0;
  prolongation[0](12,7) = -5.0/256.0;
  prolongation[0](12,8) = 15.0/256.0;
  prolongation[0](12,9) = -5.0/256.0;
  prolongation[0](12,10) = 75.0/256.0;
  prolongation[0](12,11) = -25.0/256.0;
  prolongation[0](12,12) = 225.0/256.0;
  prolongation[0](12,13) = -75.0/256.0;
  prolongation[0](12,14) = 25.0/256.0;
  prolongation[0](12,15) = -75.0/256.0;
  prolongation[0](13,4) = 5.0/16.0;
  prolongation[0](13,8) = 1.0/16.0;
  prolongation[0](13,12) = 15.0/16.0;
  prolongation[0](13,15) = -5.0/16.0;
  prolongation[0](14,12) = 1.0;
  prolongation[0](15,6) = 1.0/16.0;
  prolongation[0](15,10) = 5.0/16.0;
  prolongation[0](15,12) = 15.0/16.0;
  prolongation[0](15,13) = -5.0/16.0;
  prolongation[1](0,0) = -1.0/16.0;
  prolongation[1](0,1) = -1.0/16.0;
  prolongation[1](0,4) = 9.0/16.0;
  prolongation[1](0,5) = 9.0/16.0;
  prolongation[1](1,1) = 1.0;
  prolongation[1](2,1) = -1.0/16.0;
  prolongation[1](2,2) = -1.0/16.0;
  prolongation[1](2,6) = 9.0/16.0;
  prolongation[1](2,7) = 9.0/16.0;
  prolongation[1](3,0) = 1.0/256.0;
  prolongation[1](3,1) = 1.0/256.0;
  prolongation[1](3,2) = 1.0/256.0;
  prolongation[1](3,3) = 1.0/256.0;
  prolongation[1](3,4) = -9.0/256.0;
  prolongation[1](3,5) = -9.0/256.0;
  prolongation[1](3,6) = -9.0/256.0;
  prolongation[1](3,7) = -9.0/256.0;
  prolongation[1](3,8) = -9.0/256.0;
  prolongation[1](3,9) = -9.0/256.0;
  prolongation[1](3,10) = -9.0/256.0;
  prolongation[1](3,11) = -9.0/256.0;
  prolongation[1](3,12) = 81.0/256.0;
  prolongation[1](3,13) = 81.0/256.0;
  prolongation[1](3,14) = 81.0/256.0;
  prolongation[1](3,15) = 81.0/256.0;
  prolongation[1](4,5) = 1.0;
  prolongation[1](5,0) = 1.0/16.0;
  prolongation[1](5,1) = 5.0/16.0;
  prolongation[1](5,4) = -5.0/16.0;
  prolongation[1](5,5) = 15.0/16.0;
  prolongation[1](6,1) = 5.0/16.0;
  prolongation[1](6,2) = 1.0/16.0;
  prolongation[1](6,6) = 15.0/16.0;
  prolongation[1](6,7) = -5.0/16.0;
  prolongation[1](7,6) = 1.0;
  prolongation[1](8,5) = -1.0/16.0;
  prolongation[1](8,9) = -1.0/16.0;
  prolongation[1](8,13) = 9.0/16.0;
  prolongation[1](8,14) = 9.0/16.0;
  prolongation[1](9,0) = -1.0/256.0;
  prolongation[1](9,1) = -5.0/256.0;
  prolongation[1](9,2) = -5.0/256.0;
  prolongation[1](9,3) = -1.0/256.0;
  prolongation[1](9,4) = 5.0/256.0;
  prolongation[1](9,5) = -15.0/256.0;
  prolongation[1](9,6) = 45.0/256.0;
  prolongation[1](9,7) = 45.0/256.0;
  prolongation[1](9,8) = 5.0/256.0;
  prolongation[1](9,9) = -15.0/256.0;
  prolongation[1](9,10) = 9.0/256.0;
  prolongation[1](9,11) = 9.0/256.0;
  prolongation[1](9,12) = -45.0/256.0;
  prolongation[1](9,13) = 135.0/256.0;
  prolongation[1](9,14) = 135.0/256.0;
  prolongation[1](9,15) = -45.0/256.0;
  prolongation[1](10,0) = -5.0/256.0;
  prolongation[1](10,1) = -5.0/256.0;
  prolongation[1](10,2) = -1.0/256.0;
  prolongation[1](10,3) = -1.0/256.0;
  prolongation[1](10,4) = 45.0/256.0;
  prolongation[1](10,5) = 45.0/256.0;
  prolongation[1](10,6) = -15.0/256.0;
  prolongation[1](10,7) = 5.0/256.0;
  prolongation[1](10,8) = 9.0/256.0;
  prolongation[1](10,9) = 9.0/256.0;
  prolongation[1](10,10) = -15.0/256.0;
  prolongation[1](10,11) = 5.0/256.0;
  prolongation[1](10,12) = 135.0/256.0;
  prolongation[1](10,13) = 135.0/256.0;
  prolongation[1](10,14) = -45.0/256.0;
  prolongation[1](10,15) = -45.0/256.0;
  prolongation[1](11,6) = -1.0/16.0;
  prolongation[1](11,10) = -1.0/16.0;
  prolongation[1](11,12) = 9.0/16.0;
  prolongation[1](11,13) = 9.0/16.0;
  prolongation[1](12,5) = 5.0/16.0;
  prolongation[1](12,9) = 1.0/16.0;
  prolongation[1](12,13) = 15.0/16.0;
  prolongation[1](12,14) = -5.0/16.0;
  prolongation[1](13,0) = 5.0/256.0;
  prolongation[1](13,1) = 25.0/256.0;
  prolongation[1](13,2) = 5.0/256.0;
  prolongation[1](13,3) = 1.0/256.0;
  prolongation[1](13,4) = -25.0/256.0;
  prolongation[1](13,5) = 75.0/256.0;
  prolongation[1](13,6) = 75.0/256.0;
  prolongation[1](13,7) = -25.0/256.0;
  prolongation[1](13,8) = -5.0/256.0;
  prolongation[1](13,9) = 15.0/256.0;
  prolongation[1](13,10) = 15.0/256.0;
  prolongation[1](13,11) = -5.0/256.0;
  prolongation[1](13,12) = -75.0/256.0;
  prolongation[1](13,13) = 225.0/256.0;
  prolongation[1](13,14) = -75.0/256.0;
  prolongation[1](13,15) = 25.0/256.0;
  prolongation[1](14,6) = 5.0/16.0;
  prolongation[1](14,10) = 1.0/16.0;
  prolongation[1](14,12) = -5.0/16.0;
  prolongation[1](14,13) = 15.0/16.0;
  prolongation[1](15,13) = 1.0;
  prolongation[2](0,0) = 1.0/256.0;
  prolongation[2](0,1) = 1.0/256.0;
  prolongation[2](0,2) = 1.0/256.0;
  prolongation[2](0,3) = 1.0/256.0;
  prolongation[2](0,4) = -9.0/256.0;
  prolongation[2](0,5) = -9.0/256.0;
  prolongation[2](0,6) = -9.0/256.0;
  prolongation[2](0,7) = -9.0/256.0;
  prolongation[2](0,8) = -9.0/256.0;
  prolongation[2](0,9) = -9.0/256.0;
  prolongation[2](0,10) = -9.0/256.0;
  prolongation[2](0,11) = -9.0/256.0;
  prolongation[2](0,12) = 81.0/256.0;
  prolongation[2](0,13) = 81.0/256.0;
  prolongation[2](0,14) = 81.0/256.0;
  prolongation[2](0,15) = 81.0/256.0;
  prolongation[2](1,1) = -1.0/16.0;
  prolongation[2](1,2) = -1.0/16.0;
  prolongation[2](1,6) = 9.0/16.0;
  prolongation[2](1,7) = 9.0/16.0;
  prolongation[2](2,2) = 1.0;
  prolongation[2](3,2) = -1.0/16.0;
  prolongation[2](3,3) = -1.0/16.0;
  prolongation[2](3,8) = 9.0/16.0;
  prolongation[2](3,9) = 9.0/16.0;
  prolongation[2](4,5) = -1.0/16.0;
  prolongation[2](4,9) = -1.0/16.0;
  prolongation[2](4,13) = 9.0/16.0;
  prolongation[2](4,14) = 9.0/16.0;
  prolongation[2](5,0) = -1.0/256.0;
  prolongation[2](5,1) = -5.0/256.0;
  prolongation[2](5,2) = -5.0/256.0;
  prolongation[2](5,3) = -1.0/256.0;
  prolongation[2](5,4) = 5.0/256.0;
  prolongation[2](5,5) = -15.0/256.0;
  prolongation[2](5,6) = 45.0/256.0;
  prolongation[2](5,7) = 45.0/256.0;
  prolongation[2](5,8) = 5.0/256.0;
  prolongation[2](5,9) = -15.0/256.0;
  prolongation[2](5,10) = 9.0/256.0;
  prolongation[2](5,11) = 9.0/256.0;
  prolongation[2](5,12) = -45.0/256.0;
  prolongation[2](5,13) = 135.0/256.0;
  prolongation[2](5,14) = 135.0/256.0;
  prolongation[2](5,15) = -45.0/256.0;
  prolongation[2](6,7) = 1.0;
  prolongation[2](7,1) = 1.0/16.0;
  prolongation[2](7,2) = 5.0/16.0;
  prolongation[2](7,6) = -5.0/16.0;
  prolongation[2](7,7) = 15.0/16.0;
  prolongation[2](8,9) = 1.0;
  prolongation[2](9,2) = 5.0/16.0;
  prolongation[2](9,3) = 1.0/16.0;
  prolongation[2](9,8) = -5.0/16.0;
  prolongation[2](9,9) = 15.0/16.0;
  prolongation[2](10,7) = -1.0/16.0;
  prolongation[2](10,11) = -1.0/16.0;
  prolongation[2](10,14) = 9.0/16.0;
  prolongation[2](10,15) = 9.0/16.0;
  prolongation[2](11,0) = -1.0/256.0;
  prolongation[2](11,1) = -1.0/256.0;
  prolongation[2](11,2) = -5.0/256.0;
  prolongation[2](11,3) = -5.0/256.0;
  prolongation[2](11,4) = 9.0/256.0;
  prolongation[2](11,5) = 9.0/256.0;
  prolongation[2](11,6) = 5.0/256.0;
  prolongation[2](11,7) = -15.0/256.0;
  prolongation[2](11,8) = 45.0/256.0;
  prolongation[2](11,9) = 45.0/256.0;
  prolongation[2](11,10) = 5.0/256.0;
  prolongation[2](11,11) = -15.0/256.0;
  prolongation[2](11,12) = -45.0/256.0;
  prolongation[2](11,13) = -45.0/256.0;
  prolongation[2](11,14) = 135.0/256.0;
  prolongation[2](11,15) = 135.0/256.0;
  prolongation[2](12,14) = 1.0;
  prolongation[2](13,7) = 5.0/16.0;
  prolongation[2](13,11) = 1.0/16.0;
  prolongation[2](13,14) = 15.0/16.0;
  prolongation[2](13,15) = -5.0/16.0;
  prolongation[2](14,0) = 1.0/256.0;
  prolongation[2](14,1) = 5.0/256.0;
  prolongation[2](14,2) = 25.0/256.0;
  prolongation[2](14,3) = 5.0/256.0;
  prolongation[2](14,4) = -5.0/256.0;
  prolongation[2](14,5) = 15.0/256.0;
  prolongation[2](14,6) = -25.0/256.0;
  prolongation[2](14,7) = 75.0/256.0;
  prolongation[2](14,8) = -25.0/256.0;
  prolongation[2](14,9) = 75.0/256.0;
  prolongation[2](14,10) = -5.0/256.0;
  prolongation[2](14,11) = 15.0/256.0;
  prolongation[2](14,12) = 25.0/256.0;
  prolongation[2](14,13) = -75.0/256.0;
  prolongation[2](14,14) = 225.0/256.0;
  prolongation[2](14,15) = -75.0/256.0;
  prolongation[2](15,5) = 1.0/16.0;
  prolongation[2](15,9) = 5.0/16.0;
  prolongation[2](15,13) = -5.0/16.0;
  prolongation[2](15,14) = 15.0/16.0;
  prolongation[3](0,0) = -1.0/16.0;
  prolongation[3](0,3) = -1.0/16.0;
  prolongation[3](0,10) = 9.0/16.0;
  prolongation[3](0,11) = 9.0/16.0;
  prolongation[3](1,0) = 1.0/256.0;
  prolongation[3](1,1) = 1.0/256.0;
  prolongation[3](1,2) = 1.0/256.0;
  prolongation[3](1,3) = 1.0/256.0;
  prolongation[3](1,4) = -9.0/256.0;
  prolongation[3](1,5) = -9.0/256.0;
  prolongation[3](1,6) = -9.0/256.0;
  prolongation[3](1,7) = -9.0/256.0;
  prolongation[3](1,8) = -9.0/256.0;
  prolongation[3](1,9) = -9.0/256.0;
  prolongation[3](1,10) = -9.0/256.0;
  prolongation[3](1,11) = -9.0/256.0;
  prolongation[3](1,12) = 81.0/256.0;
  prolongation[3](1,13) = 81.0/256.0;
  prolongation[3](1,14) = 81.0/256.0;
  prolongation[3](1,15) = 81.0/256.0;
  prolongation[3](2,2) = -1.0/16.0;
  prolongation[3](2,3) = -1.0/16.0;
  prolongation[3](2,8) = 9.0/16.0;
  prolongation[3](2,9) = 9.0/16.0;
  prolongation[3](3,3) = 1.0;
  prolongation[3](4,0) = -5.0/256.0;
  prolongation[3](4,1) = -1.0/256.0;
  prolongation[3](4,2) = -1.0/256.0;
  prolongation[3](4,3) = -5.0/256.0;
  prolongation[3](4,4) = -15.0/256.0;
  prolongation[3](4,5) = 5.0/256.0;
  prolongation[3](4,6) = 9.0/256.0;
  prolongation[3](4,7) = 9.0/256.0;
  prolongation[3](4,8) = -15.0/256.0;
  prolongation[3](4,9) = 5.0/256.0;
  prolongation[3](4,10) = 45.0/256.0;
  prolongation[3](4,11) = 45.0/256.0;
  prolongation[3](4,12) = 135.0/256.0;
  prolongation[3](4,13) = -45.0/256.0;
  prolongation[3](4,14) = -45.0/256.0;
  prolongation[3](4,15) = 135.0/256.0;
  prolongation[3](5,4) = -1.0/16.0;
  prolongation[3](5,8) = -1.0/16.0;
  prolongation[3](5,12) = 9.0/16.0;
  prolongation[3](5,15) = 9.0/16.0;
  prolongation[3](6,7) = -1.0/16.0;
  prolongation[3](6,11) = -1.0/16.0;
  prolongation[3](6,14) = 9.0/16.0;
  prolongation[3](6,15) = 9.0/16.0;
  prolongation[3](7,0) = -1.0/256.0;
  prolongation[3](7,1) = -1.0/256.0;
  prolongation[3](7,2) = -5.0/256.0;
  prolongation[3](7,3) = -5.0/256.0;
  prolongation[3](7,4) = 9.0/256.0;
  prolongation[3](7,5) = 9.0/256.0;
  prolongation[3](7,6) = 5.0/256.0;
  prolongation[3](7,7) = -15.0/256.0;
  prolongation[3](7,8) = 45.0/256.0;
  prolongation[3](7,9) = 45.0/256.0;
  prolongation[3](7,10) = 5.0/256.0;
  prolongation[3](7,11) = -15.0/256.0;
  prolongation[3](7,12) = -45.0/256.0;
  prolongation[3](7,13) = -45.0/256.0;
  prolongation[3](7,14) = 135.0/256.0;
  prolongation[3](7,15) = 135.0/256.0;
  prolongation[3](8,2) = 1.0/16.0;
  prolongation[3](8,3) = 5.0/16.0;
  prolongation[3](8,8) = 15.0/16.0;
  prolongation[3](8,9) = -5.0/16.0;
  prolongation[3](9,8) = 1.0;
  prolongation[3](10,11) = 1.0;
  prolongation[3](11,0) = 1.0/16.0;
  prolongation[3](11,3) = 5.0/16.0;
  prolongation[3](11,10) = -5.0/16.0;
  prolongation[3](11,11) = 15.0/16.0;
  prolongation[3](12,7) = 1.0/16.0;
  prolongation[3](12,11) = 5.0/16.0;
  prolongation[3](12,14) = -5.0/16.0;
  prolongation[3](12,15) = 15.0/16.0;
  prolongation[3](13,15) = 1.0;
  prolongation[3](14,4) = 1.0/16.0;
  prolongation[3](14,8) = 5.0/16.0;
  prolongation[3](14,12) = -5.0/16.0;
  prolongation[3](14,15) = 15.0/16.0;
  prolongation[3](15,0) = 5.0/256.0;
  prolongation[3](15,1) = 1.0/256.0;
  prolongation[3](15,2) = 5.0/256.0;
  prolongation[3](15,3) = 25.0/256.0;
  prolongation[3](15,4) = 15.0/256.0;
  prolongation[3](15,5) = -5.0/256.0;
  prolongation[3](15,6) = -5.0/256.0;
  prolongation[3](15,7) = 15.0/256.0;
  prolongation[3](15,8) = 75.0/256.0;
  prolongation[3](15,9) = -25.0/256.0;
  prolongation[3](15,10) = -25.0/256.0;
  prolongation[3](15,11) = 75.0/256.0;
  prolongation[3](15,12) = -75.0/256.0;
  prolongation[3](15,13) = 25.0/256.0;
  prolongation[3](15,14) = -75.0/256.0;
  prolongation[3](15,15) = 225.0/256.0;

  restriction[0](0,0) = 1.0;
  restriction[0](4,5) = 1.0;
  restriction[0](10,11) = 1.0;
  restriction[0](12,14) = 1.0;
  restriction[1](1,1) = 1.0;
  restriction[1](5,4) = 1.0;
  restriction[1](6,7) = 1.0;
  restriction[1](13,15) = 1.0;
  restriction[2](2,2) = 1.0;
  restriction[2](7,6) = 1.0;
  restriction[2](9,8) = 1.0;
  restriction[2](14,12) = 1.0;
  restriction[3](3,3) = 1.0;
  restriction[3](8,9) = 1.0;
  restriction[3](11,10) = 1.0;
  restriction[3](15,13) = 1.0;
};


template <>
double
FECubicSub<2>::shape_value (const unsigned int i,
				const Point<2>    &p) const
{
  Assert (i<total_dofs, ExcInvalidIndex(i));

  const double xi = p(0),
	       eta= p(1);
  switch (i)
    {
      case 0: return 1.0-11.0/2.0*xi+9.0*xi*xi-9.0/2.0*xi*xi*xi+(-11.0/2.0+
121.0/4.0*xi-99.0/2.0*xi*xi+99.0/4.0*xi*xi*xi)*eta+(9.0-99.0/2.0*xi+81.0*xi*xi
-81.0/2.0*xi*xi*xi)*eta*eta+(-9.0/2.0+99.0/4.0*xi-81.0/2.0*xi*xi+81.0/4.0*xi*xi
*xi)*eta*eta*eta;
      case 1: return xi-9.0/2.0*xi*xi+9.0/2.0*xi*xi*xi+(-11.0/2.0*xi+99.0/4.0
*xi*xi-99.0/4.0*xi*xi*xi)*eta+(9.0*xi-81.0/2.0*xi*xi+81.0/2.0*xi*xi*xi)*eta*eta
+(-9.0/2.0*xi+81.0/4.0*xi*xi-81.0/4.0*xi*xi*xi)*eta*eta*eta;
      case 2: return (xi-9.0/2.0*xi*xi+9.0/2.0*xi*xi*xi)*eta+(-9.0/2.0*xi+
81.0/4.0*xi*xi-81.0/4.0*xi*xi*xi)*eta*eta+(9.0/2.0*xi-81.0/4.0*xi*xi+81.0/4.0*
xi*xi*xi)*eta*eta*eta;
      case 3: return (1.0-11.0/2.0*xi+9.0*xi*xi-9.0/2.0*xi*xi*xi)*eta+(-9.0/
2.0+99.0/4.0*xi-81.0/2.0*xi*xi+81.0/4.0*xi*xi*xi)*eta*eta+(9.0/2.0-99.0/4.0*xi+
81.0/2.0*xi*xi-81.0/4.0*xi*xi*xi)*eta*eta*eta;
      case 4: return 9.0*xi-45.0/2.0*xi*xi+27.0/2.0*xi*xi*xi+(-99.0/2.0*xi+
495.0/4.0*xi*xi-297.0/4.0*xi*xi*xi)*eta+(81.0*xi-405.0/2.0*xi*xi+243.0/2.0*xi*
xi*xi)*eta*eta+(-81.0/2.0*xi+405.0/4.0*xi*xi-243.0/4.0*xi*xi*xi)*eta*eta*eta;
      case 5: return -9.0/2.0*xi+18.0*xi*xi-27.0/2.0*xi*xi*xi+(99.0/4.0*xi
-99.0*xi*xi+297.0/4.0*xi*xi*xi)*eta+(-81.0/2.0*xi+162.0*xi*xi-243.0/2.0*xi*xi*
xi)*eta*eta+(81.0/4.0*xi-81.0*xi*xi+243.0/4.0*xi*xi*xi)*eta*eta*eta;
      case 6: return (9.0*xi-81.0/2.0*xi*xi+81.0/2.0*xi*xi*xi)*eta+(-45.0/2.0
*xi+405.0/4.0*xi*xi-405.0/4.0*xi*xi*xi)*eta*eta+(27.0/2.0*xi-243.0/4.0*xi*xi+
243.0/4.0*xi*xi*xi)*eta*eta*eta;
      case 7: return (-9.0/2.0*xi+81.0/4.0*xi*xi-81.0/4.0*xi*xi*xi)*eta+(18.0
*xi-81.0*xi*xi+81.0*xi*xi*xi)*eta*eta+(-27.0/2.0*xi+243.0/4.0*xi*xi-243.0/4.0*
xi*xi*xi)*eta*eta*eta;
      case 8: return (9.0*xi-45.0/2.0*xi*xi+27.0/2.0*xi*xi*xi)*eta+(-81.0/2.0
*xi+405.0/4.0*xi*xi-243.0/4.0*xi*xi*xi)*eta*eta+(81.0/2.0*xi-405.0/4.0*xi*xi+
243.0/4.0*xi*xi*xi)*eta*eta*eta;
      case 9: return (-9.0/2.0*xi+18.0*xi*xi-27.0/2.0*xi*xi*xi)*eta+(81.0/4.0
*xi-81.0*xi*xi+243.0/4.0*xi*xi*xi)*eta*eta+(-81.0/4.0*xi+81.0*xi*xi-243.0/4.0*
xi*xi*xi)*eta*eta*eta;
      case 10: return (9.0-99.0/2.0*xi+81.0*xi*xi-81.0/2.0*xi*xi*xi)*eta+(
-45.0/2.0+495.0/4.0*xi-405.0/2.0*xi*xi+405.0/4.0*xi*xi*xi)*eta*eta+(27.0/2.0
-297.0/4.0*xi+243.0/2.0*xi*xi-243.0/4.0*xi*xi*xi)*eta*eta*eta;
      case 11: return (-9.0/2.0+99.0/4.0*xi-81.0/2.0*xi*xi+81.0/4.0*xi*xi*xi)
*eta+(18.0-99.0*xi+162.0*xi*xi-81.0*xi*xi*xi)*eta*eta+(-27.0/2.0+297.0/4.0*xi
-243.0/2.0*xi*xi+243.0/4.0*xi*xi*xi)*eta*eta*eta;
      case 12: return (81.0*xi-405.0/2.0*xi*xi+243.0/2.0*xi*xi*xi)*eta+(
-405.0/2.0*xi+2025.0/4.0*xi*xi-1215.0/4.0*xi*xi*xi)*eta*eta+(243.0/2.0*xi
-1215.0/4.0*xi*xi+729.0/4.0*xi*xi*xi)*eta*eta*eta;
      case 13: return (-81.0/2.0*xi+162.0*xi*xi-243.0/2.0*xi*xi*xi)*eta+(
405.0/4.0*xi-405.0*xi*xi+1215.0/4.0*xi*xi*xi)*eta*eta+(-243.0/4.0*xi+243.0*xi*
xi-729.0/4.0*xi*xi*xi)*eta*eta*eta;
      case 14: return (81.0/4.0*xi-81.0*xi*xi+243.0/4.0*xi*xi*xi)*eta+(-81.0*
xi+324.0*xi*xi-243.0*xi*xi*xi)*eta*eta+(243.0/4.0*xi-243.0*xi*xi+729.0/4.0*xi*
xi*xi)*eta*eta*eta;
      case 15: return (-81.0/2.0*xi+405.0/4.0*xi*xi-243.0/4.0*xi*xi*xi)*eta+(
162.0*xi-405.0*xi*xi+243.0*xi*xi*xi)*eta*eta+(-243.0/2.0*xi+1215.0/4.0*xi*xi
-729.0/4.0*xi*xi*xi)*eta*eta*eta;
    };
  return 0;
};



template <>
Point<2>
FECubicSub<2>::shape_grad (const unsigned int i,
			       const Point<2>    &p) const
{
  Assert (i<total_dofs, ExcInvalidIndex(i));

  const double xi = p(0),
	       eta= p(1);
  switch (i)
    {
      case 0: return Point<2>(-11.0/2.0+18.0*xi-27.0/2.0*xi*xi+(121.0/4.0-99.0*xi+297.0/4.0*xi*xi)*eta+(-99.0/2.0+162.0*xi-243.0/2.0*xi*xi)*eta*eta+(99.0/4.0-81.0*xi+243.0/4.0*xi*xi)*eta*eta*eta,
      -11.0/2.0+121.0/4.0*xi-99.0/2.0*xi*xi+99.0/4.0*xi*xi*xi+2.0*(9.0-99.0/2.0*xi+81.0*xi*xi-81.0/2.0*xi*xi*xi)*eta+3.0*(-9.0/2.0+99.0/4.0*xi-81.0/2.0*xi*xi+81.0/4.0*xi*xi*xi)*eta*eta);
      case 1: return Point<2>(1.0-9.0*xi+27.0/2.0*xi*xi+(-11.0/2.0+99.0/2.0*xi-297.0/4.0*xi*xi)*eta+(9.0-81.0*xi+243.0/2.0*xi*xi)*eta*eta+(-9.0/2.0+81.0/2.0*xi-243.0/4.0*xi*xi)*eta*eta*eta,
      -11.0/2.0*xi+99.0/4.0*xi*xi-99.0/4.0*xi*xi*xi+2.0*(9.0*xi-81.0/2.0*xi*xi+81.0/2.0*xi*xi*xi)*eta+3.0*(-9.0/2.0*xi+81.0/4.0*xi*xi-81.0/4.0*xi*xi*xi)*eta*eta);
      case 2: return Point<2>((1.0-9.0*xi+27.0/2.0*xi*xi)*eta+(-9.0/2.0+81.0/2.0*xi-243.0/4.0*xi*xi)*eta*eta+(9.0/2.0-81.0/2.0*xi+243.0/4.0*xi*xi)*eta*eta*eta,
      xi-9.0/2.0*xi*xi+9.0/2.0*xi*xi*xi+2.0*(-9.0/2.0*xi+81.0/4.0*xi*xi-81.0/4.0*xi*xi*xi)*eta+3.0*(9.0/2.0*xi-81.0/4.0*xi*xi+81.0/4.0*xi*xi*xi)*eta*eta);
      case 3: return Point<2>((-11.0/2.0+18.0*xi-27.0/2.0*xi*xi)*eta+(99.0/4.0-81.0*xi+243.0/4.0*xi*xi)*eta*eta+(-99.0/4.0+81.0*xi-243.0/4.0*xi*xi)*eta*eta*eta,
      1.0-11.0/2.0*xi+9.0*xi*xi-9.0/2.0*xi*xi*xi+2.0*(-9.0/2.0+99.0/4.0*xi-81.0/2.0*xi*xi+81.0/4.0*xi*xi*xi)*eta+3.0*(9.0/2.0-99.0/4.0*xi+81.0/2.0*xi*xi-81.0/4.0*xi*xi*xi)*eta*eta);
      case 4: return Point<2>(9.0-45.0*xi+81.0/2.0*xi*xi+(-99.0/2.0+495.0/2.0*xi-891.0/4.0*xi*xi)*eta+(81.0-405.0*xi+729.0/2.0*xi*xi)*eta*eta+(-81.0/2.0+405.0/2.0*xi-729.0/4.0*xi*xi)*eta*eta*eta,
      -99.0/2.0*xi+495.0/4.0*xi*xi-297.0/4.0*xi*xi*xi+2.0*(81.0*xi-405.0/2.0*xi*xi+243.0/2.0*xi*xi*xi)*eta+3.0*(-81.0/2.0*xi+405.0/4.0*xi*xi-243.0/4.0*xi*xi*xi)*eta*eta);
      case 5: return Point<2>(-9.0/2.0+36.0*xi-81.0/2.0*xi*xi+(99.0/4.0-198.0*xi+891.0/4.0*xi*xi)*eta+(-81.0/2.0+324.0*xi-729.0/2.0*xi*xi)*eta*eta+(81.0/4.0-162.0*xi+729.0/4.0*xi*xi)*eta*eta*eta,
      99.0/4.0*xi-99.0*xi*xi+297.0/4.0*xi*xi*xi+2.0*(-81.0/2.0*xi+162.0*xi*xi-243.0/2.0*xi*xi*xi)*eta+3.0*(81.0/4.0*xi-81.0*xi*xi+243.0/4.0*xi*xi*xi)*eta*eta);
      case 6: return Point<2>((9.0-81.0*xi+243.0/2.0*xi*xi)*eta+(-45.0/2.0+405.0/2.0*xi-1215.0/4.0*xi*xi)*eta*eta+(27.0/2.0-243.0/2.0*xi+729.0/4.0*xi*xi)*eta*eta*eta,
      9.0*xi-81.0/2.0*xi*xi+81.0/2.0*xi*xi*xi+2.0*(-45.0/2.0*xi+405.0/4.0*xi*xi-405.0/4.0*xi*xi*xi)*eta+3.0*(27.0/2.0*xi-243.0/4.0*xi*xi+243.0/4.0*xi*xi*xi)*eta*eta);
      case 7: return Point<2>((-9.0/2.0+81.0/2.0*xi-243.0/4.0*xi*xi)*eta+(18.0-162.0*xi+243.0*xi*xi)*eta*eta+(-27.0/2.0+243.0/2.0*xi-729.0/4.0*xi*xi)*eta*eta*eta,
      -9.0/2.0*xi+81.0/4.0*xi*xi-81.0/4.0*xi*xi*xi+2.0*(18.0*xi-81.0*xi*xi+81.0*xi*xi*xi)*eta+3.0*(-27.0/2.0*xi+243.0/4.0*xi*xi-243.0/4.0*xi*xi*xi)*eta*eta);
      case 8: return Point<2>((9.0-45.0*xi+81.0/2.0*xi*xi)*eta+(-81.0/2.0+405.0/2.0*xi-729.0/4.0*xi*xi)*eta*eta+(81.0/2.0-405.0/2.0*xi+729.0/4.0*xi*xi)*eta*eta*eta,
      9.0*xi-45.0/2.0*xi*xi+27.0/2.0*xi*xi*xi+2.0*(-81.0/2.0*xi+405.0/4.0*xi*xi-243.0/4.0*xi*xi*xi)*eta+3.0*(81.0/2.0*xi-405.0/4.0*xi*xi+243.0/4.0*xi*xi*xi)*eta*eta);
      case 9: return Point<2>((-9.0/2.0+36.0*xi-81.0/2.0*xi*xi)*eta+(81.0/4.0-162.0*xi+729.0/4.0*xi*xi)*eta*eta+(-81.0/4.0+162.0*xi-729.0/4.0*xi*xi)*eta*eta*eta,
      -9.0/2.0*xi+18.0*xi*xi-27.0/2.0*xi*xi*xi+2.0*(81.0/4.0*xi-81.0*xi*xi+243.0/4.0*xi*xi*xi)*eta+3.0*(-81.0/4.0*xi+81.0*xi*xi-243.0/4.0*xi*xi*xi)*eta*eta);
      case 10: return Point<2>((-99.0/2.0+162.0*xi-243.0/2.0*xi*xi)*eta+(495.0/4.0-405.0*xi+1215.0/4.0*xi*xi)*eta*eta+(-297.0/4.0+243.0*xi-729.0/4.0*xi*xi)*eta*eta*eta,
      9.0-99.0/2.0*xi+81.0*xi*xi-81.0/2.0*xi*xi*xi+2.0*(-45.0/2.0+495.0/4.0*xi-405.0/2.0*xi*xi+405.0/4.0*xi*xi*xi)*eta+3.0*(27.0/2.0-297.0/4.0*xi+243.0/2.0*xi*xi-243.0/4.0*xi*xi*xi)*eta*eta);
      case 11: return Point<2>((99.0/4.0-81.0*xi+243.0/4.0*xi*xi)*eta+(-99.0+324.0*xi-243.0*xi*xi)*eta*eta+(297.0/4.0-243.0*xi+729.0/4.0*xi*xi)*eta*eta*eta,
      -9.0/2.0+99.0/4.0*xi-81.0/2.0*xi*xi+81.0/4.0*xi*xi*xi+2.0*(18.0-99.0*xi+162.0*xi*xi-81.0*xi*xi*xi)*eta+3.0*(-27.0/2.0+297.0/4.0*xi-243.0/2.0*xi*xi+243.0/4.0*xi*xi*xi)*eta*eta);
      case 12: return Point<2>((81.0-405.0*xi+729.0/2.0*xi*xi)*eta+(-405.0/2.0+2025.0/2.0*xi-3645.0/4.0*xi*xi)*eta*eta+(243.0/2.0-1215.0/2.0*xi+2187.0/4.0*xi*xi)*eta*eta*eta,
      81.0*xi-405.0/2.0*xi*xi+243.0/2.0*xi*xi*xi+2.0*(-405.0/2.0*xi+2025.0/4.0*xi*xi-1215.0/4.0*xi*xi*xi)*eta+3.0*(243.0/2.0*xi-1215.0/4.0*xi*xi+729.0/4.0*xi*xi*xi)*eta*eta);
      case 13: return Point<2>((-81.0/2.0+324.0*xi-729.0/2.0*xi*xi)*eta+(405.0/4.0-810.0*xi+3645.0/4.0*xi*xi)*eta*eta+(-243.0/4.0+486.0*xi-2187.0/4.0*xi*xi)*eta*eta*eta,
      -81.0/2.0*xi+162.0*xi*xi-243.0/2.0*xi*xi*xi+2.0*(405.0/4.0*xi-405.0*xi*xi+1215.0/4.0*xi*xi*xi)*eta+3.0*(-243.0/4.0*xi+243.0*xi*xi-729.0/4.0*xi*xi*xi)*eta*eta);
      case 14: return Point<2>((81.0/4.0-162.0*xi+729.0/4.0*xi*xi)*eta+(-81.0+648.0*xi-729.0*xi*xi)*eta*eta+(243.0/4.0-486.0*xi+2187.0/4.0*xi*xi)*eta*eta*eta,
      81.0/4.0*xi-81.0*xi*xi+243.0/4.0*xi*xi*xi+2.0*(-81.0*xi+324.0*xi*xi-243.0*xi*xi*xi)*eta+3.0*(243.0/4.0*xi-243.0*xi*xi+729.0/4.0*xi*xi*xi)*eta*eta);
      case 15: return Point<2>((-81.0/2.0+405.0/2.0*xi-729.0/4.0*xi*xi)*eta+(162.0-810.0*xi+729.0*xi*xi)*eta*eta+(-243.0/2.0+1215.0/2.0*xi-2187.0/4.0*xi*xi)*eta*eta*eta,
      -81.0/2.0*xi+405.0/4.0*xi*xi-243.0/4.0*xi*xi*xi+2.0*(162.0*xi-405.0*xi*xi+243.0*xi*xi*xi)*eta+3.0*(-243.0/2.0*xi+1215.0/4.0*xi*xi-729.0/4.0*xi*xi*xi)*eta*eta);
    };
  return Point<2> ();
};



template <>
Tensor<2,2>
FECubicSub<2>::shape_grad_grad (const unsigned int i,
				const Point<2>    &p) const
{
  Assert (i<total_dofs, ExcInvalidIndex(i));

  const double xi = p(0),
	       eta= p(1);
  Tensor<2,2> return_value;
  
  switch (i)
    {
      case 0:
	    return_value[0][0] = 18.0-27.0*xi+(-99.0+297.0/2.0*xi)*eta+(162.0-243.0*xi)*eta*eta+(-81.0+243.0/2.0*xi)*eta*eta*eta;
	    return_value[0][1] = 121.0/4.0-99.0*xi+297.0/4.0*xi*xi+2.0*(-99.0/2.0+162.0*xi-243.0/2.0*xi*xi)*eta+3.0*(99.0/4.0-81.0*xi+243.0/4.0*xi*xi)*eta*eta;
	    return_value[1][0] = 121.0/4.0-99.0*xi+297.0/4.0*xi*xi+2.0*(-99.0/2.0+162.0*xi-243.0/2.0*xi*xi)*eta+3.0*(99.0/4.0-81.0*xi+243.0/4.0*xi*xi)*eta*eta;
	    return_value[1][1] = 18.0-99.0*xi+162.0*xi*xi-81.0*xi*xi*xi+6.0*(-9.0/2.0+99.0/4.0*xi-81.0/2.0*xi*xi+81.0/4.0*xi*xi*xi)*eta;
	    break;
      case 1:
	    return_value[0][0] = -9.0+27.0*xi+(99.0/2.0-297.0/2.0*xi)*eta+(-81.0+243.0*xi)*eta*eta+(81.0/2.0-243.0/2.0*xi)*eta*eta*eta;
	    return_value[0][1] = -11.0/2.0+99.0/2.0*xi-297.0/4.0*xi*xi+2.0*(9.0-81.0*xi+243.0/2.0*xi*xi)*eta+3.0*(-9.0/2.0+81.0/2.0*xi-243.0/4.0*xi*xi)*eta*eta;
	    return_value[1][0] = -11.0/2.0+99.0/2.0*xi-297.0/4.0*xi*xi+2.0*(9.0-81.0*xi+243.0/2.0*xi*xi)*eta+3.0*(-9.0/2.0+81.0/2.0*xi-243.0/4.0*xi*xi)*eta*eta;
	    return_value[1][1] = 18.0*xi-81.0*xi*xi+81.0*xi*xi*xi+6.0*(-9.0/2.0*xi+81.0/4.0*xi*xi-81.0/4.0*xi*xi*xi)*eta;
	    break;
      case 2:
	    return_value[0][0] = (-9.0+27.0*xi)*eta+(81.0/2.0-243.0/2.0*xi)*eta*eta+(-81.0/2.0+243.0/2.0*xi)*eta*eta*eta;
	    return_value[0][1] = 1.0-9.0*xi+27.0/2.0*xi*xi+2.0*(-9.0/2.0+81.0/2.0*xi-243.0/4.0*xi*xi)*eta+3.0*(9.0/2.0-81.0/2.0*xi+243.0/4.0*xi*xi)*eta*eta;
	    return_value[1][0] = 1.0-9.0*xi+27.0/2.0*xi*xi+2.0*(-9.0/2.0+81.0/2.0*xi-243.0/4.0*xi*xi)*eta+3.0*(9.0/2.0-81.0/2.0*xi+243.0/4.0*xi*xi)*eta*eta;
	    return_value[1][1] = -9.0*xi+81.0/2.0*xi*xi-81.0/2.0*xi*xi*xi+6.0*(9.0/2.0*xi-81.0/4.0*xi*xi+81.0/4.0*xi*xi*xi)*eta;
	    break;
      case 3:
	    return_value[0][0] = (18.0-27.0*xi)*eta+(-81.0+243.0/2.0*xi)*eta*eta+(81.0-243.0/2.0*xi)*eta*eta*eta;
	    return_value[0][1] = -11.0/2.0+18.0*xi-27.0/2.0*xi*xi+2.0*(99.0/4.0-81.0*xi+243.0/4.0*xi*xi)*eta+3.0*(-99.0/4.0+81.0*xi-243.0/4.0*xi*xi)*eta*eta;
	    return_value[1][0] = -11.0/2.0+18.0*xi-27.0/2.0*xi*xi+2.0*(99.0/4.0-81.0*xi+243.0/4.0*xi*xi)*eta+3.0*(-99.0/4.0+81.0*xi-243.0/4.0*xi*xi)*eta*eta;
	    return_value[1][1] = -9.0+99.0/2.0*xi-81.0*xi*xi+81.0/2.0*xi*xi*xi+6.0*(9.0/2.0-99.0/4.0*xi+81.0/2.0*xi*xi-81.0/4.0*xi*xi*xi)*eta;
	    break;
      case 4:
	    return_value[0][0] = -45.0+81.0*xi+(495.0/2.0-891.0/2.0*xi)*eta+(-405.0+729.0*xi)*eta*eta+(405.0/2.0-729.0/2.0*xi)*eta*eta*eta;
	    return_value[0][1] = -99.0/2.0+495.0/2.0*xi-891.0/4.0*xi*xi+2.0*(81.0-405.0*xi+729.0/2.0*xi*xi)*eta+3.0*(-81.0/2.0+405.0/2.0*xi-729.0/4.0*xi*xi)*eta*eta;
	    return_value[1][0] = -99.0/2.0+495.0/2.0*xi-891.0/4.0*xi*xi+2.0*(81.0-405.0*xi+729.0/2.0*xi*xi)*eta+3.0*(-81.0/2.0+405.0/2.0*xi-729.0/4.0*xi*xi)*eta*eta;
	    return_value[1][1] = 162.0*xi-405.0*xi*xi+243.0*xi*xi*xi+6.0*(-81.0/2.0*xi+405.0/4.0*xi*xi-243.0/4.0*xi*xi*xi)*eta;
	    break;
      case 5:
	    return_value[0][0] = 36.0-81.0*xi+(-198.0+891.0/2.0*xi)*eta+(324.0-729.0*xi)*eta*eta+(-162.0+729.0/2.0*xi)*eta*eta*eta;
	    return_value[0][1] = 99.0/4.0-198.0*xi+891.0/4.0*xi*xi+2.0*(-81.0/2.0+324.0*xi-729.0/2.0*xi*xi)*eta+3.0*(81.0/4.0-162.0*xi+729.0/4.0*xi*xi)*eta*eta;
	    return_value[1][0] = 99.0/4.0-198.0*xi+891.0/4.0*xi*xi+2.0*(-81.0/2.0+324.0*xi-729.0/2.0*xi*xi)*eta+3.0*(81.0/4.0-162.0*xi+729.0/4.0*xi*xi)*eta*eta;
	    return_value[1][1] = -81.0*xi+324.0*xi*xi-243.0*xi*xi*xi+6.0*(81.0/4.0*xi-81.0*xi*xi+243.0/4.0*xi*xi*xi)*eta;
	    break;
      case 6:
	    return_value[0][0] = (-81.0+243.0*xi)*eta+(405.0/2.0-1215.0/2.0*xi)*eta*eta+(-243.0/2.0+729.0/2.0*xi)*eta*eta*eta;
	    return_value[0][1] = 9.0-81.0*xi+243.0/2.0*xi*xi+2.0*(-45.0/2.0+405.0/2.0*xi-1215.0/4.0*xi*xi)*eta+3.0*(27.0/2.0-243.0/2.0*xi+729.0/4.0*xi*xi)*eta*eta;
	    return_value[1][0] = 9.0-81.0*xi+243.0/2.0*xi*xi+2.0*(-45.0/2.0+405.0/2.0*xi-1215.0/4.0*xi*xi)*eta+3.0*(27.0/2.0-243.0/2.0*xi+729.0/4.0*xi*xi)*eta*eta;
	    return_value[1][1] = -45.0*xi+405.0/2.0*xi*xi-405.0/2.0*xi*xi*xi+6.0*(27.0/2.0*xi-243.0/4.0*xi*xi+243.0/4.0*xi*xi*xi)*eta;
	    break;
      case 7:
	    return_value[0][0] = (81.0/2.0-243.0/2.0*xi)*eta+(-162.0+486.0*xi)*eta*eta+(243.0/2.0-729.0/2.0*xi)*eta*eta*eta;
	    return_value[0][1] = -9.0/2.0+81.0/2.0*xi-243.0/4.0*xi*xi+2.0*(18.0-162.0*xi+243.0*xi*xi)*eta+3.0*(-27.0/2.0+243.0/2.0*xi-729.0/4.0*xi*xi)*eta*eta;
	    return_value[1][0] = -9.0/2.0+81.0/2.0*xi-243.0/4.0*xi*xi+2.0*(18.0-162.0*xi+243.0*xi*xi)*eta+3.0*(-27.0/2.0+243.0/2.0*xi-729.0/4.0*xi*xi)*eta*eta;
	    return_value[1][1] = 36.0*xi-162.0*xi*xi+162.0*xi*xi*xi+6.0*(-27.0/2.0*xi+243.0/4.0*xi*xi-243.0/4.0*xi*xi*xi)*eta;
	    break;
      case 8:
	    return_value[0][0] = (-45.0+81.0*xi)*eta+(405.0/2.0-729.0/2.0*xi)*eta*eta+(-405.0/2.0+729.0/2.0*xi)*eta*eta*eta;
	    return_value[0][1] = 9.0-45.0*xi+81.0/2.0*xi*xi+2.0*(-81.0/2.0+405.0/2.0*xi-729.0/4.0*xi*xi)*eta+3.0*(81.0/2.0-405.0/2.0*xi+729.0/4.0*xi*xi)*eta*eta;
	    return_value[1][0] = 9.0-45.0*xi+81.0/2.0*xi*xi+2.0*(-81.0/2.0+405.0/2.0*xi-729.0/4.0*xi*xi)*eta+3.0*(81.0/2.0-405.0/2.0*xi+729.0/4.0*xi*xi)*eta*eta;
	    return_value[1][1] = -81.0*xi+405.0/2.0*xi*xi-243.0/2.0*xi*xi*xi+6.0*(81.0/2.0*xi-405.0/4.0*xi*xi+243.0/4.0*xi*xi*xi)*eta;
	    break;
      case 9:
	    return_value[0][0] = (36.0-81.0*xi)*eta+(-162.0+729.0/2.0*xi)*eta*eta+(162.0-729.0/2.0*xi)*eta*eta*eta;
	    return_value[0][1] = -9.0/2.0+36.0*xi-81.0/2.0*xi*xi+2.0*(81.0/4.0-162.0*xi+729.0/4.0*xi*xi)*eta+3.0*(-81.0/4.0+162.0*xi-729.0/4.0*xi*xi)*eta*eta;
	    return_value[1][0] = -9.0/2.0+36.0*xi-81.0/2.0*xi*xi+2.0*(81.0/4.0-162.0*xi+729.0/4.0*xi*xi)*eta+3.0*(-81.0/4.0+162.0*xi-729.0/4.0*xi*xi)*eta*eta;
	    return_value[1][1] = 81.0/2.0*xi-162.0*xi*xi+243.0/2.0*xi*xi*xi+6.0*(-81.0/4.0*xi+81.0*xi*xi-243.0/4.0*xi*xi*xi)*eta;
	    return_value[0][0] = (162.0-243.0*xi)*eta+(-405.0+1215.0/2.0*xi)*eta*eta+(243.0-729.0/2.0*xi)*eta*eta*eta;
	    return_value[0][1] = -99.0/2.0+162.0*xi-243.0/2.0*xi*xi+2.0*(495.0/4.0-405.0*xi+1215.0/4.0*xi*xi)*eta+3.0*(-297.0/4.0+243.0*xi-729.0/4.0*xi*xi)*eta*eta;
	    return_value[1][0] = -99.0/2.0+162.0*xi-243.0/2.0*xi*xi+2.0*(495.0/4.0-405.0*xi+1215.0/4.0*xi*xi)*eta+3.0*(-297.0/4.0+243.0*xi-729.0/4.0*xi*xi)*eta*eta;
	    return_value[1][1] = -45.0+495.0/2.0*xi-405.0*xi*xi+405.0/2.0*xi*xi*xi+6.0*(27.0/2.0-297.0/4.0*xi+243.0/2.0*xi*xi-243.0/4.0*xi*xi*xi)*eta;
	    return_value[0][0] = (-81.0+243.0/2.0*xi)*eta+(324.0-486.0*xi)*eta*eta+(-243.0+729.0/2.0*xi)*eta*eta*eta;
	    return_value[0][1] = 99.0/4.0-81.0*xi+243.0/4.0*xi*xi+2.0*(-99.0+324.0*xi-243.0*xi*xi)*eta+3.0*(297.0/4.0-243.0*xi+729.0/4.0*xi*xi)*eta*eta;
	    return_value[1][0] = 99.0/4.0-81.0*xi+243.0/4.0*xi*xi+2.0*(-99.0+324.0*xi-243.0*xi*xi)*eta+3.0*(297.0/4.0-243.0*xi+729.0/4.0*xi*xi)*eta*eta;
	    return_value[1][1] = 36.0-198.0*xi+324.0*xi*xi-162.0*xi*xi*xi+6.0*(-27.0/2.0+297.0/4.0*xi-243.0/2.0*xi*xi+243.0/4.0*xi*xi*xi)*eta;
	    return_value[0][0] = (-405.0+729.0*xi)*eta+(2025.0/2.0-3645.0/2.0*xi)*eta*eta+(-1215.0/2.0+2187.0/2.0*xi)*eta*eta*eta;
	    return_value[0][1] = 81.0-405.0*xi+729.0/2.0*xi*xi+2.0*(-405.0/2.0+2025.0/2.0*xi-3645.0/4.0*xi*xi)*eta+3.0*(243.0/2.0-1215.0/2.0*xi+2187.0/4.0*xi*xi)*eta*eta;
	    return_value[1][0] = 81.0-405.0*xi+729.0/2.0*xi*xi+2.0*(-405.0/2.0+2025.0/2.0*xi-3645.0/4.0*xi*xi)*eta+3.0*(243.0/2.0-1215.0/2.0*xi+2187.0/4.0*xi*xi)*eta*eta;
	    return_value[1][1] = -405.0*xi+2025.0/2.0*xi*xi-1215.0/2.0*xi*xi*xi+6.0*(243.0/2.0*xi-1215.0/4.0*xi*xi+729.0/4.0*xi*xi*xi)*eta;
	    return_value[0][0] = (324.0-729.0*xi)*eta+(-810.0+3645.0/2.0*xi)*eta*eta+(486.0-2187.0/2.0*xi)*eta*eta*eta;
	    return_value[0][1] = -81.0/2.0+324.0*xi-729.0/2.0*xi*xi+2.0*(405.0/4.0-810.0*xi+3645.0/4.0*xi*xi)*eta+3.0*(-243.0/4.0+486.0*xi-2187.0/4.0*xi*xi)*eta*eta;
	    return_value[1][0] = -81.0/2.0+324.0*xi-729.0/2.0*xi*xi+2.0*(405.0/4.0-810.0*xi+3645.0/4.0*xi*xi)*eta+3.0*(-243.0/4.0+486.0*xi-2187.0/4.0*xi*xi)*eta*eta;
	    return_value[1][1] = 405.0/2.0*xi-810.0*xi*xi+1215.0/2.0*xi*xi*xi+6.0*(-243.0/4.0*xi+243.0*xi*xi-729.0/4.0*xi*xi*xi)*eta;
	    return_value[0][0] = (-162.0+729.0/2.0*xi)*eta+(648.0-1458.0*xi)*eta*eta+(-486.0+2187.0/2.0*xi)*eta*eta*eta;
	    return_value[0][1] = 81.0/4.0-162.0*xi+729.0/4.0*xi*xi+2.0*(-81.0+648.0*xi-729.0*xi*xi)*eta+3.0*(243.0/4.0-486.0*xi+2187.0/4.0*xi*xi)*eta*eta;
	    return_value[1][0] = 81.0/4.0-162.0*xi+729.0/4.0*xi*xi+2.0*(-81.0+648.0*xi-729.0*xi*xi)*eta+3.0*(243.0/4.0-486.0*xi+2187.0/4.0*xi*xi)*eta*eta;
	    return_value[1][1] = -162.0*xi+648.0*xi*xi-486.0*xi*xi*xi+6.0*(243.0/4.0*xi-243.0*xi*xi+729.0/4.0*xi*xi*xi)*eta;
	    return_value[0][0] = (405.0/2.0-729.0/2.0*xi)*eta+(-810.0+1458.0*xi)*eta*eta+(1215.0/2.0-2187.0/2.0*xi)*eta*eta*eta;
	    return_value[0][1] = -81.0/2.0+405.0/2.0*xi-729.0/4.0*xi*xi+2.0*(162.0-810.0*xi+729.0*xi*xi)*eta+3.0*(-243.0/2.0+1215.0/2.0*xi-2187.0/4.0*xi*xi)*eta*eta;
	    return_value[1][0] = -81.0/2.0+405.0/2.0*xi-729.0/4.0*xi*xi+2.0*(162.0-810.0*xi+729.0*xi*xi)*eta+3.0*(-243.0/2.0+1215.0/2.0*xi-2187.0/4.0*xi*xi)*eta*eta;
	    return_value[1][1] = 324.0*xi-810.0*xi*xi+486.0*xi*xi*xi+6.0*(-243.0/2.0*xi+1215.0/4.0*xi*xi-729.0/4.0*xi*xi*xi)*eta;
      break;
    };
  return return_value;
};



template <>
void FECubicSub<2>::get_local_mass_matrix (const DoFHandler<2>::cell_iterator &cell,
					       const Boundary<2> &,
					       dFMatrix &local_mass_matrix) const {
  Assert (local_mass_matrix.n() == total_dofs,
	  ExcWrongFieldDimension(local_mass_matrix.n(),total_dofs));
  Assert (local_mass_matrix.m() == total_dofs,
	  ExcWrongFieldDimension(local_mass_matrix.m(),total_dofs));

  const double x[4] = { cell->vertex(0)(0),
			cell->vertex(1)(0),
			cell->vertex(2)(0),
			cell->vertex(3)(0)  };
  const double y[4] = { cell->vertex(0)(1),
			cell->vertex(1)(1),
			cell->vertex(2)(1),
			cell->vertex(3)(1)  };
  
/* check that the Jacobi determinant

    t0 = (-x[0]*(1.0-eta)+x[1]*(1.0-eta)+x[2]*eta-x[3]*eta) *
         (-y[0]*(1.0-xi)-y[1]*xi+y[2]*xi+y[3]*(1.0-xi))        -
	 (-x[0]*(1.0-xi)-x[1]*xi+x[2]*xi+x[3]*(1.0-xi)) *
	 (-y[0]*(1.0-eta)+y[1]*(1.0-eta)+y[2]*eta-y[3]*eta)

   has the right sign.	 
	 
   We do not attempt to check its (hopefully) positive sign at all points
   on the unit cell, but we check that it is positive in the four corners,
   which is sufficient since $det J$ is a bilinear function.
*/
  Assert ((-x[0]+x[1])*(-y[0]+y[3])-(-x[0]+x[3])*(-y[0]+y[1]),  // xi=eta=0
	  ExcJacobiDeterminantHasWrongSign());
  Assert ((x[2]-x[3])*(-y[0]+y[3])-(-x[0]+x[3])*(y[2]-y[3]),    // xi=0, eta=1
	  ExcJacobiDeterminantHasWrongSign());
  Assert ((x[2]-x[3])*(-y[1]+y[2])-(-x[1]+x[2])*(y[2]-y[3]),    // xi=eta=1
	  ExcJacobiDeterminantHasWrongSign());
  Assert ((-x[0]+x[1])*(-y[1]+y[2])-(-x[1]+x[2])*(-y[0]+y[1]),  // xi=1, eta=0
	  ExcJacobiDeterminantHasWrongSign());

  const double t1 = x[0]-x[1]+x[2]-x[3];
  const double t2 = -y[0]+y[1];
  const double t3 = t1*t2;
  const double t4 = 19.0/44100.0*t3;
  const double t5 = -x[0]+x[1];
  const double t6 = y[0]-y[1]+y[2]-y[3];
  const double t7 = t5*t6;
  const double t8 = 19.0/44100.0*t7;
  const double t9 = -x[0]+x[3];
  const double t10 = t9*t6;
  const double t11 = 19.0/44100.0*t10;
  const double t12 = t9*t2;
  const double t13 = 64.0/11025.0*t12;
  const double t14 = -y[0]+y[3];
  const double t15 = t1*t14;
  const double t16 = 19.0/44100.0*t15;
  const double t17 = t5*t14;
  const double t18 = 64.0/11025.0*t17;
  const double t20 = 361.0/5644800.0*t10;
  const double t21 = 19.0/22050.0*t12;
  const double t22 = 361.0/5644800.0*t15;
  const double t23 = 19.0/22050.0*t17;
  const double t24 = -t4+t8-t20-t21+t22+t23;
  const double t25 = 361.0/5644800.0*t3;
  const double t26 = 361.0/5644800.0*t7;
  const double t29 = -t25+t26-t20-361.0/2822400.0*t12+t22+361.0/2822400.0*t17;
  const double t30 = -t25+t26-t11-t21+t16+t23;
  const double t31 = t3/4900.0;
  const double t32 = t7/4900.0;
  const double t33 = 209.0/627200.0*t10;
  const double t34 = 11.0/2450.0*t12;
  const double t35 = 209.0/627200.0*t15;
  const double t36 = 11.0/2450.0*t17;
  const double t37 = -t31+t32-t33-t34+t35+t36;
  const double t38 = 19.0/156800.0*t10;
  const double t39 = 2.0/1225.0*t12;
  const double t40 = 19.0/156800.0*t15;
  const double t41 = 2.0/1225.0*t17;
  const double t42 = -t31+t32+t38+t39-t40-t41;
  const double t43 = 209.0/627200.0*t3;
  const double t44 = 209.0/627200.0*t7;
  const double t45 = 19.0/627200.0*t10;
  const double t46 = 209.0/313600.0*t12;
  const double t47 = 19.0/627200.0*t15;
  const double t48 = 209.0/313600.0*t17;
  const double t49 = -t43+t44-t45-t46+t47+t48;
  const double t50 = 19.0/156800.0*t3;
  const double t51 = 19.0/156800.0*t7;
  const double t52 = 19.0/78400.0*t12;
  const double t53 = 19.0/78400.0*t17;
  const double t54 = t50-t51-t45+t52+t47-t53;
  const double t55 = 19.0/627200.0*t3;
  const double t56 = 19.0/627200.0*t7;
  const double t57 = -t55+t56-t33-t46+t35+t48;
  const double t58 = -t55+t56+t38+t52-t40-t53;
  const double t59 = t10/4900.0;
  const double t60 = t15/4900.0;
  const double t61 = -t43+t44-t59-t34+t60+t36;
  const double t62 = t50-t51-t59+t39+t60-t41;
  const double t63 = 99.0/627200.0*t3;
  const double t64 = 99.0/627200.0*t7;
  const double t65 = 99.0/627200.0*t10;
  const double t66 = 1089.0/313600.0*t12;
  const double t67 = 99.0/627200.0*t15;
  const double t68 = 1089.0/313600.0*t17;
  const double t69 = -t63+t64-t65-t66+t67+t68;
  const double t70 = 9.0/156800.0*t10;
  const double t71 = 99.0/78400.0*t12;
  const double t72 = 9.0/156800.0*t15;
  const double t73 = 99.0/78400.0*t17;
  const double t74 = -t63+t64+t70+t71-t72-t73;
  const double t75 = 9.0/156800.0*t3;
  const double t76 = 9.0/156800.0*t7;
  const double t77 = 9.0/19600.0*t12;
  const double t78 = 9.0/19600.0*t17;
  const double t79 = t75-t76+t70-t77-t72+t78;
  const double t80 = t75-t76-t65+t71+t67-t73;
  const double t81 = 79.0/14700.0*t3;
  const double t82 = 79.0/14700.0*t7;
  const double t86 = -1501.0/1881600.0*t3+1501.0/1881600.0*t7-t11-t21+t16+t23;
  const double t87 = 9.0/4900.0*t3;
  const double t88 = 9.0/4900.0*t7;
  const double t89 = t87-t88+t38+t39-t40-t41;
  const double t90 = 3.0/700.0*t3;
  const double t91 = 3.0/700.0*t7;
  const double t92 = -t90+t91-t33-t34+t35+t36;
  const double t93 = 2607.0/627200.0*t3;
  const double t94 = 2607.0/627200.0*t7;
  const double t95 = -t93+t94-t59-t34+t60+t36;
  const double t96 = 237.0/156800.0*t3;
  const double t97 = 237.0/156800.0*t7;
  const double t98 = t96-t97-t59+t39+t60-t41;
  const double t99 = 171.0/627200.0*t3;
  const double t100 = 171.0/627200.0*t7;
  const double t101 = t99-t100+t38+t52-t40-t53;
  const double t104 = -57.0/89600.0*t3+57.0/89600.0*t7-t33-t46+t35+t48;
  const double t105 = 891.0/627200.0*t3;
  const double t106 = 891.0/627200.0*t7;
  const double t107 = t105-t106+t70+t71-t72-t73;
  const double t108 = 297.0/89600.0*t3;
  const double t109 = 297.0/89600.0*t7;
  const double t110 = -t108+t109-t65-t66+t67+t68;
  const double t111 = 27.0/22400.0*t3;
  const double t112 = 27.0/22400.0*t7;
  const double t113 = t111-t112-t65+t71+t67-t73;
  const double t114 = 81.0/156800.0*t3;
  const double t115 = 81.0/156800.0*t7;
  const double t116 = -t114+t115+t70-t77-t72+t78;
  const double t117 = 79.0/14700.0*t10;
  const double t118 = 79.0/14700.0*t15;
  const double t122 = -t4+t8-1501.0/1881600.0*t10-t21+1501.0/1881600.0*t15+t23;
  const double t123 = 9.0/4900.0*t10;
  const double t124 = 9.0/4900.0*t15;
  const double t125 = t96-t97+t123+t39-t124-t41;
  const double t126 = 3.0/700.0*t10;
  const double t127 = 3.0/700.0*t15;
  const double t128 = -t93+t94-t126-t34+t127+t36;
  const double t129 = 237.0/156800.0*t10;
  const double t130 = 237.0/156800.0*t15;
  const double t131 = t87-t88+t129+t39-t130-t41;
  const double t132 = 2607.0/627200.0*t10;
  const double t133 = 2607.0/627200.0*t15;
  const double t134 = -t90+t91-t132-t34+t133+t36;
  const double t135 = 171.0/627200.0*t10;
  const double t136 = 171.0/627200.0*t15;
  const double t137 = t50-t51+t135+t52-t136-t53;
  const double t140 = -t43+t44-57.0/89600.0*t10-t46+57.0/89600.0*t15+t48;
  const double t141 = 81.0/156800.0*t10;
  const double t142 = 81.0/156800.0*t15;
  const double t143 = -t114+t115-t141-t77+t142+t78;
  const double t144 = 891.0/627200.0*t10;
  const double t145 = 891.0/627200.0*t15;
  const double t146 = t111-t112+t144+t71-t145-t73;
  const double t147 = 297.0/89600.0*t10;
  const double t148 = 297.0/89600.0*t15;
  const double t149 = -t108+t109-t147-t66+t148+t68;
  const double t150 = 27.0/22400.0*t10;
  const double t151 = 27.0/22400.0*t15;
  const double t152 = t105-t106+t150+t71-t151-t73;
  const double t154 = -t31+t32-t132-t34+t133+t36;
  const double t155 = -t31+t32+t129+t39-t130-t41;
  const double t156 = t50-t51+t123+t39-t124-t41;
  const double t157 = -t43+t44-t126-t34+t127+t36;
  const double t158 = t75-t76+t144+t71-t145-t73;
  const double t159 = t75-t76-t141-t77+t142+t78;
  const double t160 = -t63+t64+t150+t71-t151-t73;
  const double t161 = -t63+t64-t147-t66+t148+t68;
  const double t162 = 9.0/980.0*t3;
  const double t163 = 9.0/980.0*t7;
  const double t164 = 171.0/78400.0*t10;
  const double t165 = 36.0/1225.0*t12;
  const double t166 = 171.0/78400.0*t15;
  const double t167 = 36.0/1225.0*t17;
  const double t169 = 9.0/2450.0*t12;
  const double t170 = 9.0/2450.0*t17;
  const double t171 = t87-t88+t135+t169-t136-t170;
  const double t174 = 171.0/39200.0*t12;
  const double t175 = 171.0/39200.0*t17;
  const double t176 = -171.0/125440.0*t3+171.0/125440.0*t7-t164-t174+t166+t175;
  const double t179 = t99-t100+t135+171.0/313600.0*t12-t136-171.0/313600.0*t17;
  const double t180 = 891.0/125440.0*t3;
  const double t181 = 891.0/125440.0*t7;
  const double t182 = 81.0/78400.0*t10;
  const double t183 = 891.0/39200.0*t12;
  const double t184 = 81.0/78400.0*t15;
  const double t185 = 891.0/39200.0*t17;
  const double t186 = -t180+t181-t182-t183+t184+t185;
  const double t187 = 81.0/627200.0*t10;
  const double t188 = 891.0/313600.0*t12;
  const double t189 = 81.0/627200.0*t15;
  const double t190 = 891.0/313600.0*t17;
  const double t191 = t105-t106+t187+t188-t189-t190;
  const double t192 = 81.0/78400.0*t12;
  const double t193 = 81.0/78400.0*t17;
  const double t194 = -t114+t115+t187-t192-t189+t193;
  const double t195 = 81.0/31360.0*t3;
  const double t196 = 81.0/31360.0*t7;
  const double t197 = 81.0/9800.0*t12;
  const double t198 = 81.0/9800.0*t17;
  const double t199 = t195-t196-t182+t197+t184-t198;
  const double t200 = 99.0/4900.0*t3;
  const double t201 = 99.0/4900.0*t7;
  const double t205 = -1881.0/627200.0*t3+1881.0/627200.0*t7-t164-t174+t166+t175;
  const double t206 = 9801.0/627200.0*t3;
  const double t207 = 9801.0/627200.0*t7;
  const double t208 = -t206+t207-t182-t183+t184+t185;
  const double t209 = 891.0/156800.0*t3;
  const double t210 = 891.0/156800.0*t7;
  const double t211 = t209-t210-t182+t197+t184-t198;
  const double t212 = 2133.0/78400.0*t3;
  const double t213 = 2133.0/78400.0*t7;
  const double t214 = 9.0/980.0*t10;
  const double t215 = 9.0/980.0*t15;
  const double t219 = 2133.0/627200.0*t3-2133.0/627200.0*t7+t123+t169-t124-t170;
  const double t220 = 171.0/78400.0*t3;
  const double t221 = 171.0/78400.0*t7;
  const double t224 = -t220+t221-171.0/125440.0*t10-t174+171.0/125440.0*t15+t175;
  const double t225 = 729.0/78400.0*t3;
  const double t226 = 729.0/78400.0*t7;
  const double t227 = 81.0/31360.0*t10;
  const double t228 = 81.0/31360.0*t15;
  const double t229 = t225-t226+t227+t197-t228-t198;
  const double t230 = 243.0/11200.0*t3;
  const double t231 = 243.0/11200.0*t7;
  const double t232 = 891.0/125440.0*t10;
  const double t233 = 891.0/125440.0*t15;
  const double t234 = -t230+t231-t232-t183+t233+t185;
  const double t237 = 243.0/89600.0*t3-243.0/89600.0*t7+t144+t188-t145-t190;
  const double t238 = 729.0/627200.0*t3;
  const double t239 = 729.0/627200.0*t7;
  const double t240 = -t238+t239-t141-t192+t142+t193;
  const double t241 = 99.0/4900.0*t10;
  const double t242 = 99.0/4900.0*t15;
  const double t246 = -t220+t221-1881.0/627200.0*t10-t174+1881.0/627200.0*t15+t175;
  const double t247 = 9801.0/627200.0*t10;
  const double t248 = 9801.0/627200.0*t15;
  const double t249 = -t230+t231-t247-t183+t248+t185;
  const double t250 = 891.0/156800.0*t10;
  const double t251 = 891.0/156800.0*t15;
  const double t252 = t225-t226+t250+t197-t251-t198;
  const double t253 = 2133.0/78400.0*t10;
  const double t254 = 2133.0/78400.0*t15;
  const double t258 = t87-t88+2133.0/627200.0*t10+t169-2133.0/627200.0*t15-t170;
  const double t259 = 729.0/78400.0*t10;
  const double t260 = 729.0/78400.0*t15;
  const double t261 = t195-t196+t259+t197-t260-t198;
  const double t262 = 729.0/627200.0*t10;
  const double t263 = 729.0/627200.0*t15;
  const double t264 = -t114+t115-t262-t192+t263+t193;
  const double t267 = t105-t106+243.0/89600.0*t10+t188-243.0/89600.0*t15-t190;
  const double t268 = 243.0/11200.0*t10;
  const double t269 = 243.0/11200.0*t15;
  const double t270 = -t180+t181-t268-t183+t269+t185;
  const double t272 = t209-t210+t259+t197-t260-t198;
  const double t273 = -t206+t207-t268-t183+t269+t185;
  const double t275 = t99-t100+t123+t169-t124-t170;
  const double t276 = 81.0/78400.0*t3;
  const double t277 = 81.0/78400.0*t7;
  const double t278 = -t276+t277-t232-t183+t233+t185;
  const double t279 = -t276+t277+t227+t197-t228-t198;
  const double t280 = 81.0/627200.0*t3;
  const double t281 = 81.0/627200.0*t7;
  const double t282 = t280-t281-t141-t192+t142+t193;
  const double t283 = t280-t281+t144+t188-t145-t190;
  const double t285 = -t276+t277+t250+t197-t251-t198;
  const double t286 = -t276+t277-t247-t183+t248+t185;
  const double t287 = 729.0/15680.0*t3;
  const double t288 = 729.0/15680.0*t7;
  const double t289 = 729.0/15680.0*t10;
  const double t290 = 729.0/4900.0*t12;
  const double t291 = 729.0/15680.0*t15;
  const double t292 = 729.0/4900.0*t17;
  const double t295 = 729.0/39200.0*t12;
  const double t297 = 729.0/39200.0*t17;
  const double t298 = t225-t226+729.0/125440.0*t10+t295-729.0/125440.0*t15-t297;
  const double t301 = -t238+t239-t262-729.0/313600.0*t12+t263+729.0/313600.0*t17;
  const double t304 = 729.0/125440.0*t3-729.0/125440.0*t7+t259+t295-t260-t297;
  const double t305 = 8019.0/78400.0*t3;
  const double t306 = 8019.0/78400.0*t7;
  const double t310 = 8019.0/627200.0*t3-8019.0/627200.0*t7+t259+t295-t260-t297;
  const double t311 = 8019.0/78400.0*t10;
  const double t312 = 8019.0/78400.0*t15;
  const double t316 = t225-t226+8019.0/627200.0*t10+t295-8019.0/627200.0*t15-t297;
  local_mass_matrix(0,0) = -t4+t8-t11-t13+t16+t18;
  local_mass_matrix(0,1) = t24;
  local_mass_matrix(0,2) = t29;
  local_mass_matrix(0,3) = t30;
  local_mass_matrix(0,4) = t37;
  local_mass_matrix(0,5) = t42;
  local_mass_matrix(0,6) = t49;
  local_mass_matrix(0,7) = t54;
  local_mass_matrix(0,8) = t57;
  local_mass_matrix(0,9) = t58;
  local_mass_matrix(0,10) = t61;
  local_mass_matrix(0,11) = t62;
  local_mass_matrix(0,12) = t69;
  local_mass_matrix(0,13) = t74;
  local_mass_matrix(0,14) = t79;
  local_mass_matrix(0,15) = t80;
  local_mass_matrix(1,0) = t24;
  local_mass_matrix(1,1) = -t81+t82-t11-t13+t16+t18;
  local_mass_matrix(1,2) = t86;
  local_mass_matrix(1,3) = t29;
  local_mass_matrix(1,4) = t89;
  local_mass_matrix(1,5) = t92;
  local_mass_matrix(1,6) = t95;
  local_mass_matrix(1,7) = t98;
  local_mass_matrix(1,8) = t101;
  local_mass_matrix(1,9) = t104;
  local_mass_matrix(1,10) = t49;
  local_mass_matrix(1,11) = t54;
  local_mass_matrix(1,12) = t107;
  local_mass_matrix(1,13) = t110;
  local_mass_matrix(1,14) = t113;
  local_mass_matrix(1,15) = t116;
  local_mass_matrix(2,0) = t29;
  local_mass_matrix(2,1) = t86;
  local_mass_matrix(2,2) = -t81+t82-t117-t13+t118+t18;
  local_mass_matrix(2,3) = t122;
  local_mass_matrix(2,4) = t101;
  local_mass_matrix(2,5) = t104;
  local_mass_matrix(2,6) = t125;
  local_mass_matrix(2,7) = t128;
  local_mass_matrix(2,8) = t131;
  local_mass_matrix(2,9) = t134;
  local_mass_matrix(2,10) = t137;
  local_mass_matrix(2,11) = t140;
  local_mass_matrix(2,12) = t143;
  local_mass_matrix(2,13) = t146;
  local_mass_matrix(2,14) = t149;
  local_mass_matrix(2,15) = t152;
  local_mass_matrix(3,0) = t30;
  local_mass_matrix(3,1) = t29;
  local_mass_matrix(3,2) = t122;
  local_mass_matrix(3,3) = -t4+t8-t117-t13+t118+t18;
  local_mass_matrix(3,4) = t57;
  local_mass_matrix(3,5) = t58;
  local_mass_matrix(3,6) = t137;
  local_mass_matrix(3,7) = t140;
  local_mass_matrix(3,8) = t154;
  local_mass_matrix(3,9) = t155;
  local_mass_matrix(3,10) = t156;
  local_mass_matrix(3,11) = t157;
  local_mass_matrix(3,12) = t158;
  local_mass_matrix(3,13) = t159;
  local_mass_matrix(3,14) = t160;
  local_mass_matrix(3,15) = t161;
  local_mass_matrix(4,0) = t37;
  local_mass_matrix(4,1) = t89;
  local_mass_matrix(4,2) = t101;
  local_mass_matrix(4,3) = t57;
  local_mass_matrix(4,4) = -t162+t163-t164-t165+t166+t167;
  local_mass_matrix(4,5) = t171;
  local_mass_matrix(4,6) = t107;
  local_mass_matrix(4,7) = t116;
  local_mass_matrix(4,8) = t176;
  local_mass_matrix(4,9) = t179;
  local_mass_matrix(4,10) = t69;
  local_mass_matrix(4,11) = t80;
  local_mass_matrix(4,12) = t186;
  local_mass_matrix(4,13) = t191;
  local_mass_matrix(4,14) = t194;
  local_mass_matrix(4,15) = t199;
  local_mass_matrix(5,0) = t42;
  local_mass_matrix(5,1) = t92;
  local_mass_matrix(5,2) = t104;
  local_mass_matrix(5,3) = t58;
  local_mass_matrix(5,4) = t171;
  local_mass_matrix(5,5) = -t200+t201-t164-t165+t166+t167;
  local_mass_matrix(5,6) = t110;
  local_mass_matrix(5,7) = t113;
  local_mass_matrix(5,8) = t179;
  local_mass_matrix(5,9) = t205;
  local_mass_matrix(5,10) = t74;
  local_mass_matrix(5,11) = t79;
  local_mass_matrix(5,12) = t191;
  local_mass_matrix(5,13) = t208;
  local_mass_matrix(5,14) = t211;
  local_mass_matrix(5,15) = t194;
  local_mass_matrix(6,0) = t49;
  local_mass_matrix(6,1) = t95;
  local_mass_matrix(6,2) = t125;
  local_mass_matrix(6,3) = t137;
  local_mass_matrix(6,4) = t107;
  local_mass_matrix(6,5) = t110;
  local_mass_matrix(6,6) = -t212+t213-t214-t165+t215+t167;
  local_mass_matrix(6,7) = t219;
  local_mass_matrix(6,8) = t143;
  local_mass_matrix(6,9) = t146;
  local_mass_matrix(6,10) = t224;
  local_mass_matrix(6,11) = t179;
  local_mass_matrix(6,12) = t229;
  local_mass_matrix(6,13) = t234;
  local_mass_matrix(6,14) = t237;
  local_mass_matrix(6,15) = t240;
  local_mass_matrix(7,0) = t54;
  local_mass_matrix(7,1) = t98;
  local_mass_matrix(7,2) = t128;
  local_mass_matrix(7,3) = t140;
  local_mass_matrix(7,4) = t116;
  local_mass_matrix(7,5) = t113;
  local_mass_matrix(7,6) = t219;
  local_mass_matrix(7,7) = -t212+t213-t241-t165+t242+t167;
  local_mass_matrix(7,8) = t152;
  local_mass_matrix(7,9) = t149;
  local_mass_matrix(7,10) = t179;
  local_mass_matrix(7,11) = t246;
  local_mass_matrix(7,12) = t240;
  local_mass_matrix(7,13) = t237;
  local_mass_matrix(7,14) = t249;
  local_mass_matrix(7,15) = t252;
  local_mass_matrix(8,0) = t57;
  local_mass_matrix(8,1) = t101;
  local_mass_matrix(8,2) = t131;
  local_mass_matrix(8,3) = t154;
  local_mass_matrix(8,4) = t176;
  local_mass_matrix(8,5) = t179;
  local_mass_matrix(8,6) = t143;
  local_mass_matrix(8,7) = t152;
  local_mass_matrix(8,8) = -t162+t163-t253-t165+t254+t167;
  local_mass_matrix(8,9) = t258;
  local_mass_matrix(8,10) = t158;
  local_mass_matrix(8,11) = t161;
  local_mass_matrix(8,12) = t261;
  local_mass_matrix(8,13) = t264;
  local_mass_matrix(8,14) = t267;
  local_mass_matrix(8,15) = t270;
  local_mass_matrix(9,0) = t58;
  local_mass_matrix(9,1) = t104;
  local_mass_matrix(9,2) = t134;
  local_mass_matrix(9,3) = t155;
  local_mass_matrix(9,4) = t179;
  local_mass_matrix(9,5) = t205;
  local_mass_matrix(9,6) = t146;
  local_mass_matrix(9,7) = t149;
  local_mass_matrix(9,8) = t258;
  local_mass_matrix(9,9) = -t200+t201-t253-t165+t254+t167;
  local_mass_matrix(9,10) = t159;
  local_mass_matrix(9,11) = t160;
  local_mass_matrix(9,12) = t264;
  local_mass_matrix(9,13) = t272;
  local_mass_matrix(9,14) = t273;
  local_mass_matrix(9,15) = t267;
  local_mass_matrix(10,0) = t61;
  local_mass_matrix(10,1) = t49;
  local_mass_matrix(10,2) = t137;
  local_mass_matrix(10,3) = t156;
  local_mass_matrix(10,4) = t69;
  local_mass_matrix(10,5) = t74;
  local_mass_matrix(10,6) = t224;
  local_mass_matrix(10,7) = t179;
  local_mass_matrix(10,8) = t158;
  local_mass_matrix(10,9) = t159;
  local_mass_matrix(10,10) = -t220+t221-t214-t165+t215+t167;
  local_mass_matrix(10,11) = t275;
  local_mass_matrix(10,12) = t278;
  local_mass_matrix(10,13) = t279;
  local_mass_matrix(10,14) = t282;
  local_mass_matrix(10,15) = t283;
  local_mass_matrix(11,0) = t62;
  local_mass_matrix(11,1) = t54;
  local_mass_matrix(11,2) = t140;
  local_mass_matrix(11,3) = t157;
  local_mass_matrix(11,4) = t80;
  local_mass_matrix(11,5) = t79;
  local_mass_matrix(11,6) = t179;
  local_mass_matrix(11,7) = t246;
  local_mass_matrix(11,8) = t161;
  local_mass_matrix(11,9) = t160;
  local_mass_matrix(11,10) = t275;
  local_mass_matrix(11,11) = -t220+t221-t241-t165+t242+t167;
  local_mass_matrix(11,12) = t283;
  local_mass_matrix(11,13) = t282;
  local_mass_matrix(11,14) = t285;
  local_mass_matrix(11,15) = t286;
  local_mass_matrix(12,0) = t69;
  local_mass_matrix(12,1) = t107;
  local_mass_matrix(12,2) = t143;
  local_mass_matrix(12,3) = t158;
  local_mass_matrix(12,4) = t186;
  local_mass_matrix(12,5) = t191;
  local_mass_matrix(12,6) = t229;
  local_mass_matrix(12,7) = t240;
  local_mass_matrix(12,8) = t261;
  local_mass_matrix(12,9) = t264;
  local_mass_matrix(12,10) = t278;
  local_mass_matrix(12,11) = t283;
  local_mass_matrix(12,12) = -t287+t288-t289-t290+t291+t292;
  local_mass_matrix(12,13) = t298;
  local_mass_matrix(12,14) = t301;
  local_mass_matrix(12,15) = t304;
  local_mass_matrix(13,0) = t74;
  local_mass_matrix(13,1) = t110;
  local_mass_matrix(13,2) = t146;
  local_mass_matrix(13,3) = t159;
  local_mass_matrix(13,4) = t191;
  local_mass_matrix(13,5) = t208;
  local_mass_matrix(13,6) = t234;
  local_mass_matrix(13,7) = t237;
  local_mass_matrix(13,8) = t264;
  local_mass_matrix(13,9) = t272;
  local_mass_matrix(13,10) = t279;
  local_mass_matrix(13,11) = t282;
  local_mass_matrix(13,12) = t298;
  local_mass_matrix(13,13) = -t305+t306-t289-t290+t291+t292;
  local_mass_matrix(13,14) = t310;
  local_mass_matrix(13,15) = t301;
  local_mass_matrix(14,0) = t79;
  local_mass_matrix(14,1) = t113;
  local_mass_matrix(14,2) = t149;
  local_mass_matrix(14,3) = t160;
  local_mass_matrix(14,4) = t194;
  local_mass_matrix(14,5) = t211;
  local_mass_matrix(14,6) = t237;
  local_mass_matrix(14,7) = t249;
  local_mass_matrix(14,8) = t267;
  local_mass_matrix(14,9) = t273;
  local_mass_matrix(14,10) = t282;
  local_mass_matrix(14,11) = t285;
  local_mass_matrix(14,12) = t301;
  local_mass_matrix(14,13) = t310;
  local_mass_matrix(14,14) = -t305+t306-t311-t290+t312+t292;
  local_mass_matrix(14,15) = t316;
  local_mass_matrix(15,0) = t80;
  local_mass_matrix(15,1) = t116;
  local_mass_matrix(15,2) = t152;
  local_mass_matrix(15,3) = t161;
  local_mass_matrix(15,4) = t199;
  local_mass_matrix(15,5) = t194;
  local_mass_matrix(15,6) = t240;
  local_mass_matrix(15,7) = t252;
  local_mass_matrix(15,8) = t270;
  local_mass_matrix(15,9) = t267;
  local_mass_matrix(15,10) = t283;
  local_mass_matrix(15,11) = t286;
  local_mass_matrix(15,12) = t304;
  local_mass_matrix(15,13) = t301;
  local_mass_matrix(15,14) = t316;
  local_mass_matrix(15,15) = -t287+t288-t311-t290+t312+t292;
};



template <>
void FECubicSub<2>::get_unit_support_points (vector<Point<2> > &unit_points) const {
  Assert (unit_points.size() == total_dofs,
	  ExcWrongFieldDimension (unit_points.size(), total_dofs));

  unit_points[0] = Point<2>(0,0);
  unit_points[1] = Point<2>(1,0);
  unit_points[2] = Point<2>(1,1);
  unit_points[3] = Point<2>(0,1);
  unit_points[4] = Point<2>(1./3,0);
  unit_points[5] = Point<2>(2./3,0);
  unit_points[6] = Point<2>(1,1./3);
  unit_points[7] = Point<2>(1,2./3);
  unit_points[8] = Point<2>(1./3,1);
  unit_points[9] = Point<2>(2./3,1);
  unit_points[10]= Point<2>(0,1./3);
  unit_points[11]= Point<2>(0,2./3);
  unit_points[12]= Point<2>(1./3,1./3);
  unit_points[13]= Point<2>(2./3,1./3);
  unit_points[14]= Point<2>(2./3,2./3);
  unit_points[15]= Point<2>(1./3,2./3);
};



template <>
void FECubicSub<2>::get_support_points (const typename DoFHandler<2>::cell_iterator &cell,
					   const Boundary<2>&,
					   vector<Point<2> >  &support_points) const {
  Assert (support_points.size() == total_dofs,
	  ExcWrongFieldDimension (support_points.size(), total_dofs));

  const double x[4] = { cell->vertex(0)(0),
			cell->vertex(1)(0),
			cell->vertex(2)(0),
			cell->vertex(3)(0)  };
  const double y[4] = { cell->vertex(0)(1),
			cell->vertex(1)(1),
			cell->vertex(2)(1),
			cell->vertex(3)(1)  };
  const double t1 = 2.0/3.0*x[0];
  const double t2 = x[1]/3.0;
  const double t4 = 2.0/3.0*y[0];
  const double t5 = y[1]/3.0;
  const double t7 = x[0]/3.0;
  const double t8 = 2.0/3.0*x[1];
  const double t10 = y[0]/3.0;
  const double t11 = 2.0/3.0*y[1];
  const double t13 = x[2]/3.0;
  const double t15 = y[2]/3.0;
  const double t17 = 2.0/3.0*x[2];
  const double t19 = 2.0/3.0*y[2];
  const double t21 = 2.0/3.0*x[3];
  const double t23 = 2.0/3.0*y[3];
  const double t25 = x[3]/3.0;
  const double t27 = y[3]/3.0;
  const double t34 = 2.0/9.0*x[1];
  const double t36 = 2.0/9.0*x[3];
  const double t39 = 2.0/9.0*y[1];
  const double t41 = 2.0/9.0*y[3];
  const double t43 = 2.0/9.0*x[0];
  const double t45 = 2.0/9.0*x[2];
  const double t48 = 2.0/9.0*y[0];
  const double t50 = 2.0/9.0*y[2];
  support_points[0](0) = x[0];
  support_points[0](1) = y[0];
  support_points[1](0) = x[1];
  support_points[1](1) = y[1];
  support_points[2](0) = x[2];
  support_points[2](1) = y[2];
  support_points[3](0) = x[3];
  support_points[3](1) = y[3];
  support_points[4](0) = t1+t2;
  support_points[4](1) = t4+t5;
  support_points[5](0) = t7+t8;
  support_points[5](1) = t10+t11;
  support_points[6](0) = t8+t13;
  support_points[6](1) = t11+t15;
  support_points[7](0) = t2+t17;
  support_points[7](1) = t5+t19;
  support_points[8](0) = t13+t21;
  support_points[8](1) = t15+t23;
  support_points[9](0) = t17+t25;
  support_points[9](1) = t19+t27;
  support_points[10](0) = t1+t25;
  support_points[10](1) = t4+t27;
  support_points[11](0) = t7+t21;
  support_points[11](1) = t10+t23;
  support_points[12](0) = 4.0/9.0*x[0]+t34+x[2]/9.0+t36;
  support_points[12](1) = 4.0/9.0*y[0]+t39+y[2]/9.0+t41;
  support_points[13](0) = t43+4.0/9.0*x[1]+t45+x[3]/9.0;
  support_points[13](1) = t48+4.0/9.0*y[1]+t50+y[3]/9.0;
  support_points[14](0) = x[0]/9.0+t34+4.0/9.0*x[2]+t36;
  support_points[14](1) = y[0]/9.0+t39+4.0/9.0*y[2]+t41;
  support_points[15](0) = t43+x[1]/9.0+t45+4.0/9.0*x[3];
  support_points[15](1) = t48+y[1]/9.0+t50+4.0/9.0*y[3];
};



template <>
void FECubicSub<2>::get_face_support_points (const typename DoFHandler<2>::face_iterator &face,
						const Boundary<2>  &,
						vector<Point<2> >  &support_points) const {
  Assert (support_points.size() == dofs_per_face,
	  ExcWrongFieldDimension (support_points.size(), dofs_per_face));

  for (unsigned int vertex=0; vertex<2; ++vertex)
    support_points[vertex] = face->vertex(vertex);
  support_points[2] = (2*support_points[0] + support_points[1]) / 3;
  support_points[3] = (support_points[0] + 2*support_points[1]) / 3;
};



#endif






// explicit instantiations

template class FECubicSub<deal_II_dimension>;

