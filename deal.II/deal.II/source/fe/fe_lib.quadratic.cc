/* $Id$ */
/* Copyright W. Bangerth, University of Heidelberg, 1998 */


#include <fe/fe_lib.lagrange.h>
#include <grid/tria_iterator.h>
#include <grid/dof_accessor.h>
#include <grid/geometry_info.h>
#include <algorithm>



// declare explicit specializations before use:
template <> void FEQuadraticSub<deal_II_dimension>::initialize_matrices ();



#if deal_II_dimension == 1

template <>
FEQuadraticSub<1>::FEQuadraticSub () :
		FELinearMapping<1> (1, 1) {
  initialize_matrices ();
};



template <>
FEQuadraticSub<1>::FEQuadraticSub (const int) :
		FELinearMapping<1> (0, 3) {
  initialize_matrices ();
};



template <>
void FEQuadraticSub<1>::initialize_matrices () {
/*
  Get the prolongation matrices by the following little maple script:

  phi[0] := proc(xi) (1-xi)*(1-2*xi); end;
  phi[1] := proc(xi) xi*(2*xi-1);     end;
  phi[2] := proc(xi) 4*xi*(1-xi);     end;

  points[0] := array(0..2, [0, 1/2, 1/4]);
  points[1] := array(0..2, [1/2, 1, 3/4]);

  prolongation := array(0..1,0..2, 0..2);
  restriction := array(0..1,0..2, 0..2);

  for i from 0 to 1 do
    for j from 0 to 2 do
      for k from 0 to 2 do
        prolongation[i,j,k] := phi[k](points[i][j]);
      od;
    od;
  od;


  
  # to get the restriction (interpolation) matrices, evaluate
  # the basis functions on the child cells at the global
  # interpolation points

  global_points := array(0..2, [0,1,1/2]):
  child_phi[0] := proc(i, point)
                    if ((point<0) or (point>1/2)) then
		      0:
		    else
		      phi[i](2*point):
		    fi:
		  end: 
  child_phi[1] := proc(i, point)
                    if ((point<1/2) or (point>1)) then
		      0:
		    else
		      phi[i](2*point-1):
		    fi:
		  end: 
  
  for child from 0 to 1 do
    for j from 0 to 2 do
      for k from 0 to 2 do
        restriction[child,j,k] := child_phi[child](k, global_points[j]):
      od:
    od:
  od:
  
  readlib(C);
  C(prolongation);
  C(restriction);
*/

  prolongation[0](0,0) = 1.0;
  prolongation[0](0,1) = 0.0;
  prolongation[0](0,2) = 0.0;
  prolongation[0](1,0) = 0.0;
  prolongation[0](1,1) = 0.0;
  prolongation[0](1,2) = 1.0;
  prolongation[0](2,0) = 3.0/8.0;
  prolongation[0](2,1) = -1.0/8.0;
  prolongation[0](2,2) = 3.0/4.0;
  prolongation[1](0,0) = 0.0;
  prolongation[1](0,1) = 0.0;
  prolongation[1](0,2) = 1.0;
  prolongation[1](1,0) = 0.0;
  prolongation[1](1,1) = 1.0;
  prolongation[1](1,2) = 0.0;
  prolongation[1](2,0) = -1.0/8.0;
  prolongation[1](2,1) = 3.0/8.0;
  prolongation[1](2,2) = 3.0/4.0;

  restriction[0](0,0)= 1.0;
  restriction[0](2,1)= 1.0;
  restriction[1](1,1)= 1.0;
  restriction[1](2,0)= 1.0;
};



template <>
double
FEQuadraticSub<1>::shape_value(const unsigned int i,
			       const Point<1>     &p) const
{
  Assert((i<total_dofs), ExcInvalidIndex(i));
  const double xi = p(0);
  switch (i)
    {
      case 0: return (1-xi)*(1-2*xi);
      case 1: return xi*(2*xi-1);
      case 2: return 4*xi*(1-xi);
    }
  return 0.;
};



template <>
Tensor<1,1>
FEQuadraticSub<1>::shape_grad(const unsigned int i,
			      const Point<1>    &p) const
{
  Assert((i<total_dofs), ExcInvalidIndex(i));
  const double xi = p(0);
				   // originally, the return type of the
				   // function was Point<dim>, so we
				   // still construct it as that. it should
				   // make no difference in practice,
				   // however
  switch (i)
    {
      case 0: return Point<1>(-3+4*xi);
      case 1: return Point<1>(4*xi-1);
      case 2: return Point<1>(4-8*xi);
    }
  return Point<1>();
};



template <>
Tensor<2,1>
FEQuadraticSub<1>::shape_grad_grad (const unsigned int i,
				    const Point<1>    &) const
{
  Assert((i<total_dofs), ExcInvalidIndex(i));

  Tensor<2,1> return_value;
  switch (i)
    {
      case 0:
	    return_value[0][0] = 4;
	    break;
      case 1:
	    return_value[0][0] = 4;
	    break;
      case 2:
	    return_value[0][0] = -8;
	    break;
    }
  return return_value;
};



template <>
void FEQuadraticSub<1>::get_unit_support_points (vector<Point<1> > &unit_points) const {
  FiniteElement<1>::get_unit_support_points (unit_points);
};



template <>
void FEQuadraticSub<1>::get_support_points (const typename DoFHandler<1>::cell_iterator &cell,
					   const Boundary<1>  &boundary,
					   vector<Point<1> >  &support_points) const {
  FiniteElement<1>::get_support_points (cell, boundary, support_points);
};



template <>
void FEQuadraticSub<1>::get_face_support_points (const typename DoFHandler<1>::face_iterator &,
					     const Boundary<1>  &,
					     vector<Point<1> >  &) const {
  Assert (false, ExcInternalError());
};



template <>
void FEQuadraticSub<1>::get_local_mass_matrix (const DoFHandler<1>::cell_iterator &cell,
					       const Boundary<1> &,
					       dFMatrix &local_mass_matrix) const {
  Assert (local_mass_matrix.n() == total_dofs,
	  ExcWrongFieldDimension(local_mass_matrix.n(),total_dofs));
  Assert (local_mass_matrix.m() == total_dofs,
	  ExcWrongFieldDimension(local_mass_matrix.m(),total_dofs));

  const double h = cell->vertex(1)(0) - cell->vertex(0)(0);
  Assert (h>0, ExcJacobiDeterminantHasWrongSign());

  local_mass_matrix(0,0) = local_mass_matrix(1,1) = 2./15.*h;
  local_mass_matrix(0,1) = local_mass_matrix(1,0) = -1./30.*h;
  local_mass_matrix(0,2) = local_mass_matrix(2,0) = 1./15.*h;
  local_mass_matrix(1,2) = local_mass_matrix(2,1) = 1./15.*h;
  local_mass_matrix(2,2) = 8./15.*h;
};

#endif


#if deal_II_dimension == 2

template <>
FEQuadraticSub<2>::FEQuadraticSub () :
		FELinearMapping<2> (1, 1, 1)
{
  interface_constraints(0,2) = 1.0;
  interface_constraints(1,0) = 3./8.;
  interface_constraints(1,1) = -1./8.;
  interface_constraints(1,2) = 3./4.;
  interface_constraints(2,0) = -1./8.;
  interface_constraints(2,1) = 3./8.;
  interface_constraints(2,2) = 3./4.;

  initialize_matrices ();
};



template <>
FEQuadraticSub<2>::FEQuadraticSub (const int) :
		FELinearMapping<2> (0, 0, 9)
{
  initialize_matrices ();
};



template <>
void FEQuadraticSub<2>::initialize_matrices () {
/*
  Get the prolongation and restriction matrices by the following little maple script:

  phi[0] := proc(xi,eta) (1-xi)*( 2*xi-1) * (1-eta)*( 2*eta-1);    end;
  phi[1] := proc(xi,eta)    xi *(-2*xi+1) * (1-eta)*( 2*eta-1);    end;
  phi[2] := proc(xi,eta)    xi *(-2*xi+1) *    eta *(-2*eta+1);    end;
  phi[3] := proc(xi,eta) (1-xi)*( 2*xi-1) *    eta *(-2*eta+1);    end;
  phi[4] := proc(xi,eta) 4 * (1-xi)*xi        * (1-eta)*(1-2*eta); end;
  phi[5] := proc(xi,eta) 4 *    xi *(-1+2*xi) * (1-eta)*eta;       end;
  phi[6] := proc(xi,eta) 4 * (1-xi)*xi        *    eta *(-1+2*eta);end;
  phi[7] := proc(xi,eta) 4 * (1-xi)*(1-2*xi)  * (1-eta)*eta;       end;
  phi[8] := proc(xi,eta) 16 * xi*(1-xi) * eta*(1-eta);             end;

  points_x[0] := array(0..8, [0, 1/2, 1/2, 0, 1/4, 1/2, 1/4, 0, 1/4]);
  points_y[0] := array(0..8, [0, 0, 1/2, 1/2, 0, 1/4, 1/2, 1/4, 1/4]);

  points_x[1] := array(0..8, [1/2, 1, 1, 1/2, 3/4, 1, 3/4, 1/2, 3/4]);
  points_y[1] := array(0..8, [0, 0, 1/2, 1/2, 0, 1/4, 1/2, 1/4, 1/4]);

  points_x[2] := array(0..8, [1/2, 1, 1, 1/2, 3/4, 1, 3/4, 1/2, 3/4]);
  points_y[2] := array(0..8, [1/2, 1/2, 1, 1, 1/2, 3/4, 1, 3/4, 3/4]);

  points_x[3] := array(0..8, [0, 1/2, 1/2, 0, 1/4, 1/2, 1/4, 0, 1/4]);
  points_y[3] := array(0..8, [1/2, 1/2, 1, 1, 1/2, 3/4, 1, 3/4, 3/4]);

  prolongation := array(0..3,0..8, 0..8);

  for i from 0 to 3 do
    for j from 0 to 8 do
      for k from 0 to 8 do
        prolongation[i,j,k] := phi[k](points_x[i][j], points_y[i][j]);
      od;
    od;
  od;

  readlib(C);
  C(prolongation);
*/

  prolongation[0](0,0) = 1.0;
  prolongation[0](0,1) = 0.0;
  prolongation[0](0,2) = 0.0;
  prolongation[0](0,3) = 0.0;
  prolongation[0](0,4) = 0.0;
  prolongation[0](0,5) = 0.0;
  prolongation[0](0,6) = 0.0;
  prolongation[0](0,7) = 0.0;
  prolongation[0](0,8) = 0.0;
  prolongation[0](1,0) = 0.0;
  prolongation[0](1,1) = 0.0;
  prolongation[0](1,2) = 0.0;
  prolongation[0](1,3) = 0.0;
  prolongation[0](1,4) = 1.0;
  prolongation[0](1,5) = 0.0;
  prolongation[0](1,6) = 0.0;
  prolongation[0](1,7) = 0.0;
  prolongation[0](1,8) = 0.0;
  prolongation[0](2,0) = 0.0;
  prolongation[0](2,1) = 0.0;
  prolongation[0](2,2) = 0.0;
  prolongation[0](2,3) = 0.0;
  prolongation[0](2,4) = 0.0;
  prolongation[0](2,5) = 0.0;
  prolongation[0](2,6) = 0.0;
  prolongation[0](2,7) = 0.0;
  prolongation[0](2,8) = 1.0;
  prolongation[0](3,0) = 0.0;
  prolongation[0](3,1) = 0.0;
  prolongation[0](3,2) = 0.0;
  prolongation[0](3,3) = 0.0;
  prolongation[0](3,4) = 0.0;
  prolongation[0](3,5) = 0.0;
  prolongation[0](3,6) = 0.0;
  prolongation[0](3,7) = 1.0;
  prolongation[0](3,8) = 0.0;
  prolongation[0](4,0) = 3.0/8.0;
  prolongation[0](4,1) = -1.0/8.0;
  prolongation[0](4,2) = 0.0;
  prolongation[0](4,3) = 0.0;
  prolongation[0](4,4) = 3.0/4.0;
  prolongation[0](4,5) = 0.0;
  prolongation[0](4,6) = 0.0;
  prolongation[0](4,7) = 0.0;
  prolongation[0](4,8) = 0.0;
  prolongation[0](5,0) = 0.0;
  prolongation[0](5,1) = 0.0;
  prolongation[0](5,2) = 0.0;
  prolongation[0](5,3) = 0.0;
  prolongation[0](5,4) = 3.0/8.0;
  prolongation[0](5,5) = 0.0;
  prolongation[0](5,6) = -1.0/8.0;
  prolongation[0](5,7) = 0.0;
  prolongation[0](5,8) = 3.0/4.0;
  prolongation[0](6,0) = 0.0;
  prolongation[0](6,1) = 0.0;
  prolongation[0](6,2) = 0.0;
  prolongation[0](6,3) = 0.0;
  prolongation[0](6,4) = 0.0;
  prolongation[0](6,5) = -1.0/8.0;
  prolongation[0](6,6) = 0.0;
  prolongation[0](6,7) = 3.0/8.0;
  prolongation[0](6,8) = 3.0/4.0;
  prolongation[0](7,0) = 3.0/8.0;
  prolongation[0](7,1) = 0.0;
  prolongation[0](7,2) = 0.0;
  prolongation[0](7,3) = -1.0/8.0;
  prolongation[0](7,4) = 0.0;
  prolongation[0](7,5) = 0.0;
  prolongation[0](7,6) = 0.0;
  prolongation[0](7,7) = 3.0/4.0;
  prolongation[0](7,8) = 0.0;
  prolongation[0](8,0) = 9.0/64.0;
  prolongation[0](8,1) = -3.0/64.0;
  prolongation[0](8,2) = 1.0/64.0;
  prolongation[0](8,3) = -3.0/64.0;
  prolongation[0](8,4) = 9.0/32.0;
  prolongation[0](8,5) = -3.0/32.0;
  prolongation[0](8,6) = -3.0/32.0;
  prolongation[0](8,7) = 9.0/32.0;
  prolongation[0](8,8) = 9.0/16.0;
  prolongation[1](0,0) = 0.0;
  prolongation[1](0,1) = 0.0;
  prolongation[1](0,2) = 0.0;
  prolongation[1](0,3) = 0.0;
  prolongation[1](0,4) = 1.0;
  prolongation[1](0,5) = 0.0;
  prolongation[1](0,6) = 0.0;
  prolongation[1](0,7) = 0.0;
  prolongation[1](0,8) = 0.0;
  prolongation[1](1,0) = 0.0;
  prolongation[1](1,1) = 1.0;
  prolongation[1](1,2) = 0.0;
  prolongation[1](1,3) = 0.0;
  prolongation[1](1,4) = 0.0;
  prolongation[1](1,5) = 0.0;
  prolongation[1](1,6) = 0.0;
  prolongation[1](1,7) = 0.0;
  prolongation[1](1,8) = 0.0;
  prolongation[1](2,0) = 0.0;
  prolongation[1](2,1) = 0.0;
  prolongation[1](2,2) = 0.0;
  prolongation[1](2,3) = 0.0;
  prolongation[1](2,4) = 0.0;
  prolongation[1](2,5) = 1.0;
  prolongation[1](2,6) = 0.0;
  prolongation[1](2,7) = 0.0;
  prolongation[1](2,8) = 0.0;
  prolongation[1](3,0) = 0.0;
  prolongation[1](3,1) = 0.0;
  prolongation[1](3,2) = 0.0;
  prolongation[1](3,3) = 0.0;
  prolongation[1](3,4) = 0.0;
  prolongation[1](3,5) = 0.0;
  prolongation[1](3,6) = 0.0;
  prolongation[1](3,7) = 0.0;
  prolongation[1](3,8) = 1.0;
  prolongation[1](4,0) = -1.0/8.0;
  prolongation[1](4,1) = 3.0/8.0;
  prolongation[1](4,2) = 0.0;
  prolongation[1](4,3) = 0.0;
  prolongation[1](4,4) = 3.0/4.0;
  prolongation[1](4,5) = 0.0;
  prolongation[1](4,6) = 0.0;
  prolongation[1](4,7) = 0.0;
  prolongation[1](4,8) = 0.0;
  prolongation[1](5,0) = 0.0;
  prolongation[1](5,1) = 3.0/8.0;
  prolongation[1](5,2) = -1.0/8.0;
  prolongation[1](5,3) = 0.0;
  prolongation[1](5,4) = 0.0;
  prolongation[1](5,5) = 3.0/4.0;
  prolongation[1](5,6) = 0.0;
  prolongation[1](5,7) = 0.0;
  prolongation[1](5,8) = 0.0;
  prolongation[1](6,0) = 0.0;
  prolongation[1](6,1) = 0.0;
  prolongation[1](6,2) = 0.0;
  prolongation[1](6,3) = 0.0;
  prolongation[1](6,4) = 0.0;
  prolongation[1](6,5) = 3.0/8.0;
  prolongation[1](6,6) = 0.0;
  prolongation[1](6,7) = -1.0/8.0;
  prolongation[1](6,8) = 3.0/4.0;
  prolongation[1](7,0) = 0.0;
  prolongation[1](7,1) = 0.0;
  prolongation[1](7,2) = 0.0;
  prolongation[1](7,3) = 0.0;
  prolongation[1](7,4) = 3.0/8.0;
  prolongation[1](7,5) = 0.0;
  prolongation[1](7,6) = -1.0/8.0;
  prolongation[1](7,7) = 0.0;
  prolongation[1](7,8) = 3.0/4.0;
  prolongation[1](8,0) = -3.0/64.0;
  prolongation[1](8,1) = 9.0/64.0;
  prolongation[1](8,2) = -3.0/64.0;
  prolongation[1](8,3) = 1.0/64.0;
  prolongation[1](8,4) = 9.0/32.0;
  prolongation[1](8,5) = 9.0/32.0;
  prolongation[1](8,6) = -3.0/32.0;
  prolongation[1](8,7) = -3.0/32.0;
  prolongation[1](8,8) = 9.0/16.0;
  prolongation[2](0,0) = 0.0;
  prolongation[2](0,1) = 0.0;
  prolongation[2](0,2) = 0.0;
  prolongation[2](0,3) = 0.0;
  prolongation[2](0,4) = 0.0;
  prolongation[2](0,5) = 0.0;
  prolongation[2](0,6) = 0.0;
  prolongation[2](0,7) = 0.0;
  prolongation[2](0,8) = 1.0;
  prolongation[2](1,0) = 0.0;
  prolongation[2](1,1) = 0.0;
  prolongation[2](1,2) = 0.0;
  prolongation[2](1,3) = 0.0;
  prolongation[2](1,4) = 0.0;
  prolongation[2](1,5) = 1.0;
  prolongation[2](1,6) = 0.0;
  prolongation[2](1,7) = 0.0;
  prolongation[2](1,8) = 0.0;
  prolongation[2](2,0) = 0.0;
  prolongation[2](2,1) = 0.0;
  prolongation[2](2,2) = 1.0;
  prolongation[2](2,3) = 0.0;
  prolongation[2](2,4) = 0.0;
  prolongation[2](2,5) = 0.0;
  prolongation[2](2,6) = 0.0;
  prolongation[2](2,7) = 0.0;
  prolongation[2](2,8) = 0.0;
  prolongation[2](3,0) = 0.0;
  prolongation[2](3,1) = 0.0;
  prolongation[2](3,2) = 0.0;
  prolongation[2](3,3) = 0.0;
  prolongation[2](3,4) = 0.0;
  prolongation[2](3,5) = 0.0;
  prolongation[2](3,6) = 1.0;
  prolongation[2](3,7) = 0.0;
  prolongation[2](3,8) = 0.0;
  prolongation[2](4,0) = 0.0;
  prolongation[2](4,1) = 0.0;
  prolongation[2](4,2) = 0.0;
  prolongation[2](4,3) = 0.0;
  prolongation[2](4,4) = 0.0;
  prolongation[2](4,5) = 3.0/8.0;
  prolongation[2](4,6) = 0.0;
  prolongation[2](4,7) = -1.0/8.0;
  prolongation[2](4,8) = 3.0/4.0;
  prolongation[2](5,0) = 0.0;
  prolongation[2](5,1) = -1.0/8.0;
  prolongation[2](5,2) = 3.0/8.0;
  prolongation[2](5,3) = 0.0;
  prolongation[2](5,4) = 0.0;
  prolongation[2](5,5) = 3.0/4.0;
  prolongation[2](5,6) = 0.0;
  prolongation[2](5,7) = 0.0;
  prolongation[2](5,8) = 0.0;
  prolongation[2](6,0) = 0.0;
  prolongation[2](6,1) = 0.0;
  prolongation[2](6,2) = 3.0/8.0;
  prolongation[2](6,3) = -1.0/8.0;
  prolongation[2](6,4) = 0.0;
  prolongation[2](6,5) = 0.0;
  prolongation[2](6,6) = 3.0/4.0;
  prolongation[2](6,7) = 0.0;
  prolongation[2](6,8) = 0.0;
  prolongation[2](7,0) = 0.0;
  prolongation[2](7,1) = 0.0;
  prolongation[2](7,2) = 0.0;
  prolongation[2](7,3) = 0.0;
  prolongation[2](7,4) = -1.0/8.0;
  prolongation[2](7,5) = 0.0;
  prolongation[2](7,6) = 3.0/8.0;
  prolongation[2](7,7) = 0.0;
  prolongation[2](7,8) = 3.0/4.0;
  prolongation[2](8,0) = 1.0/64.0;
  prolongation[2](8,1) = -3.0/64.0;
  prolongation[2](8,2) = 9.0/64.0;
  prolongation[2](8,3) = -3.0/64.0;
  prolongation[2](8,4) = -3.0/32.0;
  prolongation[2](8,5) = 9.0/32.0;
  prolongation[2](8,6) = 9.0/32.0;
  prolongation[2](8,7) = -3.0/32.0;
  prolongation[2](8,8) = 9.0/16.0;
  prolongation[3](0,0) = 0.0;
  prolongation[3](0,1) = 0.0;
  prolongation[3](0,2) = 0.0;
  prolongation[3](0,3) = 0.0;
  prolongation[3](0,4) = 0.0;
  prolongation[3](0,5) = 0.0;
  prolongation[3](0,6) = 0.0;
  prolongation[3](0,7) = 1.0;
  prolongation[3](0,8) = 0.0;
  prolongation[3](1,0) = 0.0;
  prolongation[3](1,1) = 0.0;
  prolongation[3](1,2) = 0.0;
  prolongation[3](1,3) = 0.0;
  prolongation[3](1,4) = 0.0;
  prolongation[3](1,5) = 0.0;
  prolongation[3](1,6) = 0.0;
  prolongation[3](1,7) = 0.0;
  prolongation[3](1,8) = 1.0;
  prolongation[3](2,0) = 0.0;
  prolongation[3](2,1) = 0.0;
  prolongation[3](2,2) = 0.0;
  prolongation[3](2,3) = 0.0;
  prolongation[3](2,4) = 0.0;
  prolongation[3](2,5) = 0.0;
  prolongation[3](2,6) = 1.0;
  prolongation[3](2,7) = 0.0;
  prolongation[3](2,8) = 0.0;
  prolongation[3](3,0) = 0.0;
  prolongation[3](3,1) = 0.0;
  prolongation[3](3,2) = 0.0;
  prolongation[3](3,3) = 1.0;
  prolongation[3](3,4) = 0.0;
  prolongation[3](3,5) = 0.0;
  prolongation[3](3,6) = 0.0;
  prolongation[3](3,7) = 0.0;
  prolongation[3](3,8) = 0.0;
  prolongation[3](4,0) = 0.0;
  prolongation[3](4,1) = 0.0;
  prolongation[3](4,2) = 0.0;
  prolongation[3](4,3) = 0.0;
  prolongation[3](4,4) = 0.0;
  prolongation[3](4,5) = -1.0/8.0;
  prolongation[3](4,6) = 0.0;
  prolongation[3](4,7) = 3.0/8.0;
  prolongation[3](4,8) = 3.0/4.0;
  prolongation[3](5,0) = 0.0;
  prolongation[3](5,1) = 0.0;
  prolongation[3](5,2) = 0.0;
  prolongation[3](5,3) = 0.0;
  prolongation[3](5,4) = -1.0/8.0;
  prolongation[3](5,5) = 0.0;
  prolongation[3](5,6) = 3.0/8.0;
  prolongation[3](5,7) = 0.0;
  prolongation[3](5,8) = 3.0/4.0;
  prolongation[3](6,0) = 0.0;
  prolongation[3](6,1) = 0.0;
  prolongation[3](6,2) = -1.0/8.0;
  prolongation[3](6,3) = 3.0/8.0;
  prolongation[3](6,4) = 0.0;
  prolongation[3](6,5) = 0.0;
  prolongation[3](6,6) = 3.0/4.0;
  prolongation[3](6,7) = 0.0;
  prolongation[3](6,8) = 0.0;
  prolongation[3](7,0) = -1.0/8.0;
  prolongation[3](7,1) = 0.0;
  prolongation[3](7,2) = 0.0;
  prolongation[3](7,3) = 3.0/8.0;
  prolongation[3](7,4) = 0.0;
  prolongation[3](7,5) = 0.0;
  prolongation[3](7,6) = 0.0;
  prolongation[3](7,7) = 3.0/4.0;
  prolongation[3](7,8) = 0.0;
  prolongation[3](8,0) = -3.0/64.0;
  prolongation[3](8,1) = 1.0/64.0;
  prolongation[3](8,2) = -3.0/64.0;
  prolongation[3](8,3) = 9.0/64.0;
  prolongation[3](8,4) = -3.0/32.0;
  prolongation[3](8,5) = -3.0/32.0;
  prolongation[3](8,6) = 9.0/32.0;
  prolongation[3](8,7) = 9.0/32.0;
  prolongation[3](8,8) = 9.0/16.0;

  restriction[0](0,0) = 1.0;
  restriction[0](4,1) = 1.0;
  restriction[0](7,3) = 1.0;
  restriction[0](8,2) = 1.0;
  restriction[1](1,1) = 1.0;
  restriction[1](4,0) = 1.0;
  restriction[1](5,2) = 1.0;
  restriction[1](8,3) = 1.0;
  restriction[2](2,2) = 1.0;
  restriction[2](5,1) = 1.0;
  restriction[2](6,3) = 1.0;
  restriction[2](8,0) = 1.0;
  restriction[3](3,3) = 1.0;
  restriction[3](6,2) = 1.0;
  restriction[3](7,0) = 1.0;
  restriction[3](8,1) = 1.0; 
};



template <>
double
FEQuadraticSub<2>::shape_value (const unsigned int i,
				const Point<2>    &p) const
{
  Assert (i<total_dofs, ExcInvalidIndex(i));

  const double xi = p(0),
	       eta= p(1);
  switch (i)
    {
      case 0: return (1-xi)*( 2*xi-1) * (1-eta)*( 2*eta-1);
      case 1: return    xi *(-2*xi+1) * (1-eta)*( 2*eta-1);
      case 2: return    xi *(-2*xi+1) *    eta *(-2*eta+1);
      case 3: return (1-xi)*( 2*xi-1) *    eta *(-2*eta+1);
      case 4: return 4 * (1-xi)*xi        * (1-eta)*(1-2*eta);
      case 5: return 4 *    xi *(-1+2*xi) * (1-eta)*eta;
      case 6: return 4 * (1-xi)*xi        *    eta *(-1+2*eta);
      case 7: return 4 * (1-xi)*(1-2*xi)  * (1-eta)*eta;
      case 8: return 16 * xi*(1-xi) * eta*(1-eta);
    };
  return 0;
};



template <>
Tensor<1,2>
FEQuadraticSub<2>::shape_grad (const unsigned int i,
			       const Point<2>    &p) const
{
  Assert (i<total_dofs, ExcInvalidIndex(i));

  const double xi = p(0),
	       eta= p(1);
				   // originally, the return type of the
				   // function was Point<dim>, so we
				   // still construct it as that. it should
				   // make no difference in practice,
				   // however
  switch (i)
    {
      case 0: return Point<2>(-(2*xi-1)*(1-eta)*(2*eta-1)+2*(1-xi)*(1-eta)*(2*eta-1),
			      -(1-xi)*(2*xi-1)*(2*eta-1)+2*(1-xi)*(2*xi-1)*(1-eta));
      case 1: return Point<2>((-2*xi+1)*(1-eta)*(2*eta-1)-2*xi*(1-eta)*(2*eta-1),
			      -xi*(-2*xi+1)*(2*eta-1)+2*xi*(-2*xi+1)*(1-eta));
      case 2: return Point<2>((-2*xi+1)*eta*(-2*eta+1)-2*xi*eta*(-2*eta+1),
			      xi*(-2*xi+1)*(-2*eta+1)-2*xi*(-2*xi+1)*eta);
      case 3: return Point<2>(-(2*xi-1)*eta*(-2*eta+1)+2*(1-xi)*eta*(-2*eta+1),
			      (1-xi)*(2*xi-1)*(-2*eta+1)-2*(1-xi)*(2*xi-1)*eta);
      case 4: return Point<2>(-4*xi*(1-eta)*(-2*eta+1)+4*(1-xi)*(1-eta)*(-2*eta+1),
			      -4*(1-xi)*xi*(-2*eta+1)-8*(1-xi)*xi*(1-eta));
      case 5: return Point<2>(4*(2*xi-1)*(1-eta)*eta+8*xi*(1-eta)*eta,
			      -4*xi*(2*xi-1)*eta+4*xi*(2*xi-1)*(1-eta));
      case 6: return Point<2>(-4*xi*eta*(2*eta-1)+4*(1-xi)*eta*(2*eta-1),
			      4*(1-xi)*xi*(2*eta-1)+8*(1-xi)*xi*eta);
      case 7: return Point<2>(-4*(-2*xi+1)*(1-eta)*eta-8*(1-xi)*(1-eta)*eta,
			      -4*(1-xi)*(-2*xi+1)*eta+4*(1-xi)*(-2*xi+1)*(1-eta));
      case 8: return Point<2>(16*(1-xi)*(1-eta)*eta-16*xi*eta*(1-eta),
			      16*xi*(1-xi)*(1-eta)-16*(1-xi)*xi*eta);
    };
  return Point<2> ();
};



template <>
Tensor<2,2>
FEQuadraticSub<2>::shape_grad_grad (const unsigned int i,
				    const Point<2>    &p) const
{
  Assert (i<total_dofs, ExcInvalidIndex(i));

  const double xi = p(0),
	       eta= p(1);
  Tensor<2,2> return_value;
  
  switch (i)
    {
      case 0:
	    return_value[0][0] = 4.0-12.0*eta+8.0*eta*eta;
	    return_value[0][1] = 9.0-12.0*xi+2.0*(-6.0+8.0*xi)*eta;
	    return_value[1][0] = 9.0-12.0*xi+2.0*(-6.0+8.0*xi)*eta;
	    return_value[1][1] = 4.0-12.0*xi+8.0*xi*xi;
	    break;
      case 1:
	    return_value[0][0] = 4.0-12.0*eta+8.0*eta*eta;
	    return_value[0][1] = 3.0-12.0*xi+2.0*(-2.0+8.0*xi)*eta;
	    return_value[1][0] = 3.0-12.0*xi+2.0*(-2.0+8.0*xi)*eta;
	    return_value[1][1] = -4.0*xi+8.0*xi*xi;
	    break;
      case 2:
	    return_value[0][0] = -4.0*eta+8.0*eta*eta;
	    return_value[0][1] = 1.0-4.0*xi+2.0*(-2.0+8.0*xi)*eta;
	    return_value[1][0] = 1.0-4.0*xi+2.0*(-2.0+8.0*xi)*eta;
	    return_value[1][1] = -4.0*xi+8.0*xi*xi;
	    break;
      case 3:
	    return_value[0][0] = -4.0*eta+8.0*eta*eta;
	    return_value[0][1] = 3.0-4.0*xi+2.0*(-6.0+8.0*xi)*eta;
	    return_value[1][0] = 3.0-4.0*xi+2.0*(-6.0+8.0*xi)*eta;
	    return_value[1][1] = 4.0-12.0*xi+8.0*xi*xi;
	    break;
      case 4:
	    return_value[0][0] = -8.0+24.0*eta-16.0*eta*eta;
	    return_value[0][1] = -12.0+24.0*xi+2.0*(8.0-16.0*xi)*eta;
	    return_value[1][0] = -12.0+24.0*xi+2.0*(8.0-16.0*xi)*eta;
	    return_value[1][1] = 16.0*xi-16.0*xi*xi;
	    break;
      case 5:
	    return_value[0][0] = 16.0*eta-16.0*eta*eta;
	    return_value[0][1] = -4.0+16.0*xi+2.0*(4.0-16.0*xi)*eta;
	    return_value[1][0] = -4.0+16.0*xi+2.0*(4.0-16.0*xi)*eta;
	    return_value[1][1] = 8.0*xi-16.0*xi*xi;
	    break;
      case 6:
	    return_value[0][0] = 8.0*eta-16.0*eta*eta;
	    return_value[0][1] = -4.0+8.0*xi+2.0*(8.0-16.0*xi)*eta;
	    return_value[1][0] = -4.0+8.0*xi+2.0*(8.0-16.0*xi)*eta;
	    return_value[1][1] = 16.0*xi-16.0*xi*xi;
	    break;
      case 7:
	    return_value[0][0] = 16.0*eta-16.0*eta*eta;
	    return_value[0][1] = -12.0+16.0*xi+2.0*(12.0-16.0*xi)*eta;
	    return_value[1][0] = -12.0+16.0*xi+2.0*(12.0-16.0*xi)*eta;
	    return_value[1][1] = -8.0+24.0*xi-16.0*xi*xi;
	    break;
      case 8:
	    return_value[0][0] = -32.0*eta+32.0*eta*eta;
	    return_value[0][1] = 16.0-32.0*xi+2.0*(-16.0+32.0*xi)*eta;
	    return_value[1][0] = 16.0-32.0*xi+2.0*(-16.0+32.0*xi)*eta;
	    return_value[1][1] = -32.0*xi+32.0*xi*xi;
      break;
    };
  return return_value;
};



template <>
void FEQuadraticSub<2>::get_local_mass_matrix (const DoFHandler<2>::cell_iterator &cell,
					       const Boundary<2> &,
					       dFMatrix &local_mass_matrix) const {
  Assert (local_mass_matrix.n() == total_dofs,
	  ExcWrongFieldDimension(local_mass_matrix.n(),total_dofs));
  Assert (local_mass_matrix.m() == total_dofs,
	  ExcWrongFieldDimension(local_mass_matrix.m(),total_dofs));

/* Get the computation of the local mass matrix by these lines in maple. Note
   that tphi[i] are the basis function of the linear finite element, which
   are used by the transformation (therefore >t<phi), while the phi[i]
   are the basis functions of the biquadratic element.

   x_real := sum(x[i]*tphi[i], i=0..3);
   y_real := sum(y[i]*tphi[i], i=0..3);
   tphi[0] := (1-xi)*(1-eta);
   tphi[1] := xi*(1-eta);
   tphi[2] := xi*eta;
   tphi[3] := (1-xi)*eta;
   detJ := diff(x_real,xi)*diff(y_real,eta) - diff(x_real,eta)*diff(y_real,xi);

   phi[0] := (1-xi)*( 2*xi-1) * (1-eta)*( 2*eta-1);
   phi[1] :=    xi *(-2*xi+1) * (1-eta)*( 2*eta-1);
   phi[2] :=    xi *(-2*xi+1) *    eta *(-2*eta+1);
   phi[3] := (1-xi)*( 2*xi-1) *    eta *(-2*eta+1);
   phi[4] := 4 * (1-xi)*xi        * (1-eta)*(1-2*eta);
   phi[5] := 4 *    xi *(-1+2*xi) * (1-eta)*eta;
   phi[6] := 4 * (1-xi)*xi        *    eta *(-1+2*eta);
   phi[7] := 4 * (1-xi)*(1-2*xi)  * (1-eta)*eta;
   phi[8] := 16 * xi*(1-xi) * eta*(1-eta);
   m := proc (i,j)  int( int(phi[i]*phi[j]*detJ, xi=0..1), eta=0..1); end;

   M := array(0..8,0..8);
   for i from 0 to 8 do
     for j from 0 to 8 do
       M[i,j] := m(i,j);
     od;
   od;

   readlib(C);
   C(M, optimized);
*/

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

  const double t1 = (x[1]*y[0]);
  const double t2 = (x[1]*y[2]);
  const double t3 = (x[0]*y[3]);
  const double t4 = (x[3]*y[2]);
  const double t5 = (x[2]*y[3]);
  const double t6 = (x[0]*y[1]);
  const double t7 = (x[3]*y[1]);
  const double t8 = (x[3]*y[0]);
  const double t9 = (x[2]*y[1]);
  const double t10 = (x[1]*y[3]);
  const double t12 = (x[0]*y[2]);
  const double t13 = (x[2]*y[0]);
  const double t14 = (7.0/1800.0*t1-t2/450+t3/450+t4/1800-t5/1800-
		      7.0/1800.0*t6+t12/600+
		      t7/600-t8/450-t13/600+t9/450-t10/600);
  const double t15 = (-t1/1800+t2/1800-t3/1800-t4/1800+t5/1800+
		      t6/1800+t8/1800-t9/1800);
  const double t16 = (t1/450-t2/1800+7.0/1800.0*t3+t4/450-
		      t5/450-t6/450-t12/600+t7/600
		      -7.0/1800.0*t8+t13/600+t9/1800-t10/600);
  const double t17 = (-7.0/900.0*t1-2.0/225.0*t3-t4/900+t5/900
		      +7.0/900.0*t6+t12/900-7.0/
		      900.0*t7+2.0/225.0*t8-t13/900+7.0/900.0*t10);
  const double t18 = (t1/450-t2/900+t3/900-t6/450+t12/900+
		      t7/900-t8/900-t13/900+t9/900-
		      t10/900);
  const double t19 = (t1/900+t3/450+t4/900-t5/900-t6/900
		      -t12/900+t7/900-t8/450+t13/900-
		      t10/900);
  const double t20 = (-2.0/225.0*t1+t2/900-7.0/900.0*t3+
		      2.0/225.0*t6-t12/900-7.0/900.0*t7
		      +7.0/900.0*t8+t13/900-t9/900+7.0/900.0*t10);
  const double t21 = (-t1/225-t3/225+t6/225-t7/225+t8/225+t10/225);
  const double t23 = (t1/450-7.0/1800.0*t2+t3/1800+t4/450
		      -t5/450-t6/450+t12/600-t7/600-t8
		      /1800-t13/600+7.0/1800.0*t9+t10/600);
  const double t24 = (-7.0/900.0*t1+2.0/225.0*t2-t4/900+t5/900
		      +7.0/900.0*t6-7.0/900.0*t12
		      +t7/900+7.0/900.0*t13-2.0/225.0*t9-t10/900);
  const double t25 = (-2.0/225.0*t1+7.0/900.0*t2-t3/900+2.0/225.0*t6
		      -7.0/900.0*t12-t7/900
		      +t8/900+7.0/900.0*t13-7.0/900.0*t9+t10/900);
  const double t26 = (t1/900-t2/450+t4/900-t5/900-t6/900+t12/900
		      -t7/900-t13/900+t9/450+
		      t10/900);
  const double t27 = (-t1/225+t2/225+t6/225-t12/225+t13/225-t9/225);
  const double t29 = (t1/1800-t2/450+t3/450+7.0/1800.0*t4-7.0/1800.0*t5
		      -t6/1800-t12/600-
		      t7/600-t8/450+t13/600+t9/450+t10/600);
  const double t30 = (7.0/900.0*t2-t3/900-2.0/225.0*t4+2.0/225.0*t5
		      +t12/900+7.0/900.0*t7+
		      t8/900-t13/900-7.0/900.0*t9-7.0/900.0*t10);
  const double t31 = (-t1/900+2.0/225.0*t2-7.0/900.0*t4+7.0/900.0*t5
		      +t6/900-t12/900+7.0/
		      900.0*t7+t13/900-2.0/225.0*t9-7.0/900.0*t10);
  const double t32 = (-t2/900+t3/900+t4/450-t5/450-t12/900-t7/900
		      -t8/900+t13/900+t9/900+
		      t10/900);
  const double t33 = (t2/225-t4/225+t5/225+t7/225-t9/225-t10/225);
  const double t35 = (-t1/900-2.0/225.0*t3-7.0/900.0*t4+7.0/900.0*t5
		      +t6/900+7.0/900.0*t12
		      -t7/900+2.0/225.0*t8-7.0/900.0*t13+t10/900);
  const double t36 = (t2/900-7.0/900.0*t3-2.0/225.0*t4+2.0/225.0*t5
		      +7.0/900.0*t12+t7/900+
		      7.0/900.0*t8-7.0/900.0*t13-t9/900-t10/900);
  const double t37 = (-t3/225-t4/225+t5/225+t12/225+t8/225-t13/225);
  const double t38 = (-14.0/225.0*t1+8.0/225.0*t2-8.0/225.0*t3
		      -2.0/225.0*t4+2.0/225.0*t5+
		      14.0/225.0*t6-2.0/75.0*t12-2.0/75.0*t7
		      +8.0/225.0*t8+2.0/75.0*t13-8.0/225.0*t9+
		      2.0/75.0*t10);
  const double t39 = (2.0/225.0*t1-2.0/225.0*t2+2.0/225.0*t3
		      +2.0/225.0*t4-2.0/225.0*t5
		      -2.0/225.0*t6-2.0/225.0*t8+2.0/225.0*t9);
  const double t40 = (-8.0/225.0*t1+4.0/225.0*t2-4.0/225.0*t3
		      +8.0/225.0*t6-4.0/225.0*t12
		      -4.0/225.0*t7+4.0/225.0*t8+4.0/225.0*t13
		      -4.0/225.0*t9+4.0/225.0*t10);
  const double t41 = (-8.0/225.0*t1+14.0/225.0*t2-2.0/225.0*t3
		      -8.0/225.0*t4+8.0/225.0*t5+
		      8.0/225.0*t6-2.0/75.0*t12+2.0/75.0*t7
		      +2.0/225.0*t8+2.0/75.0*t13-14.0/225.0*t9
		      -2.0/75.0*t10);
  const double t42 = (-4.0/225.0*t1+8.0/225.0*t2-4.0/225.0*t4
		      +4.0/225.0*t5+4.0/225.0*t6
		      -4.0/225.0*t12+4.0/225.0*t7+4.0/225.0*t13
		      -8.0/225.0*t9-4.0/225.0*t10);
  const double t43 = (-2.0/225.0*t1+8.0/225.0*t2-8.0/225.0*t3
		      -14.0/225.0*t4+14.0/225.0*t5
		      +2.0/225.0*t6+2.0/75.0*t12+2.0/75.0*t7
		      +8.0/225.0*t8-2.0/75.0*t13-8.0/225.0*t9
		      -2.0/75.0*t10);
  const double t44 = (4.0/225.0*t2-4.0/225.0*t3-8.0/225.0*t4
		      +8.0/225.0*t5+4.0/225.0*t12+
		      4.0/225.0*t7+4.0/225.0*t8-4.0/225.0*t13
		      -4.0/225.0*t9-4.0/225.0*t10);
  const double t45 = (-8.0/225.0*t1+2.0/225.0*t2-14.0/225.0*t3
		      -8.0/225.0*t4+8.0/225.0*t5+
		      8.0/225.0*t6+2.0/75.0*t12-2.0/75.0*t7
		      +14.0/225.0*t8-2.0/75.0*t13-2.0/225.0*t9+
		      2.0/75.0*t10);
  const double t46 = (-4.0/225.0*t1-8.0/225.0*t3-4.0/225.0*t4
		      +4.0/225.0*t5+4.0/225.0*t6+
		      4.0/225.0*t12-4.0/225.0*t7+8.0/225.0*t8
		      -4.0/225.0*t13+4.0/225.0*t10);
  
  local_mass_matrix(0,0) = (-7.0/450.0*t1+t2/450-7.0/450.0*t3
			    -t4/450+t5/450+7.0/450.0*t6-t7/75
			    +7.0/450.0*t8-t9/450+t10/75);
  local_mass_matrix(0,1) = (t14);
  local_mass_matrix(0,2) = (t15);
  local_mass_matrix(0,3) = (t16);
  local_mass_matrix(0,4) = (t17);
  local_mass_matrix(0,5) = (t18);
  local_mass_matrix(0,6) = (t19);
  local_mass_matrix(0,7) = (t20);
  local_mass_matrix(0,8) = (t21);
  local_mass_matrix(1,0) = (t14);
  local_mass_matrix(1,1) = (-7.0/450.0*t1+7.0/450.0*t2-t3/450
			    -t4/450+t5/450+7.0/450.0*t6-
			    t12/75+t8/450+t13/75-7.0/450.0*t9);
  local_mass_matrix(1,2) = (t23);
  local_mass_matrix(1,3) = (t15);
  local_mass_matrix(1,4) = (t24);
  local_mass_matrix(1,5) = (t25);
  local_mass_matrix(1,6) = (t26);
  local_mass_matrix(1,7) = (t18);
  local_mass_matrix(1,8) = (t27);
  local_mass_matrix(2,0) = (t15);
  local_mass_matrix(2,1) = (t23);
  local_mass_matrix(2,2) = (-t1/450+7.0/450.0*t2-t3/450-7.0/450.0*t4
			    +7.0/450.0*t5+t6/450+t7/75
			    +t8/450-7.0/450.0*t9-t10/75);
  local_mass_matrix(2,3) = (t29);
  local_mass_matrix(2,4) = (t26);
  local_mass_matrix(2,5) = (t30);
  local_mass_matrix(2,6) = (t31);
  local_mass_matrix(2,7) = (t32);
  local_mass_matrix(2,8) = (t33);
  local_mass_matrix(3,0) = (t16);
  local_mass_matrix(3,1) = (t15);
  local_mass_matrix(3,2) = (t29);
  local_mass_matrix(3,3) = (-t1/450+t2/450-7.0/450.0*t3-7.0/450.0*t4
			    +7.0/450.0*t5+t6/450+
			    t12/75+7.0/450.0*t8-t13/75-t9/450);
  local_mass_matrix(3,4) = (t19);
  local_mass_matrix(3,5) = (t32);
  local_mass_matrix(3,6) = (t35);
  local_mass_matrix(3,7) = (t36);
  local_mass_matrix(3,8) = (t37);
  local_mass_matrix(4,0) = (t17);
  local_mass_matrix(4,1) = (t24);
  local_mass_matrix(4,2) = (t26);
  local_mass_matrix(4,3) = (t19);
  local_mass_matrix(4,4) = (t38);
  local_mass_matrix(4,5) = (t27);
  local_mass_matrix(4,6) = (t39);
  local_mass_matrix(4,7) = (t21);
  local_mass_matrix(4,8) = (t40);
  local_mass_matrix(5,0) = (t18);
  local_mass_matrix(5,1) = (t25);
  local_mass_matrix(5,2) = (t30);
  local_mass_matrix(5,3) = (t32);
  local_mass_matrix(5,4) = (t27);
  local_mass_matrix(5,5) = (t41);
  local_mass_matrix(5,6) = (t33);
  local_mass_matrix(5,7) = (t39);
  local_mass_matrix(5,8) = (t42);
  local_mass_matrix(6,0) = (t19);
  local_mass_matrix(6,1) = (t26);
  local_mass_matrix(6,2) = (t31);
  local_mass_matrix(6,3) = (t35);
  local_mass_matrix(6,4) = (t39);
  local_mass_matrix(6,5) = (t33);
  local_mass_matrix(6,6) = (t43);
  local_mass_matrix(6,7) = (t37);
  local_mass_matrix(6,8) = (t44);
  local_mass_matrix(7,0) = (t20);
  local_mass_matrix(7,1) = (t18);
  local_mass_matrix(7,2) = (t32);
  local_mass_matrix(7,3) = (t36);
  local_mass_matrix(7,4) = (t21);
  local_mass_matrix(7,5) = (t39);
  local_mass_matrix(7,6) = (t37);
  local_mass_matrix(7,7) = (t45);
  local_mass_matrix(7,8) = (t46);
  local_mass_matrix(8,0) = (t21);
  local_mass_matrix(8,1) = (t27);
  local_mass_matrix(8,2) = (t33);
  local_mass_matrix(8,3) = (t37);
  local_mass_matrix(8,4) = (t40);
  local_mass_matrix(8,5) = (t42);
  local_mass_matrix(8,6) = (t44);
  local_mass_matrix(8,7) = (t46);
  local_mass_matrix(8,8) = (-32.0/225.0*t1+32.0/225.0*t2-32.0/225.0*t3
			    -32.0/225.0*t4+32.0/225.0*t5+32.0/225.0*t6
			    +32.0/225.0*t8-32.0/225.0*t9);  
};



template <>
void FEQuadraticSub<2>::get_unit_support_points (vector<Point<2> > &unit_points) const {
  Assert (unit_points.size() == total_dofs,
	  ExcWrongFieldDimension (unit_points.size(), total_dofs));
  
  unit_points[0] = Point<2> (0,0);
  unit_points[1] = Point<2> (1,0);
  unit_points[2] = Point<2> (1,1);
  unit_points[3] = Point<2> (0,1);
  unit_points[4] = Point<2> (0.5,0);
  unit_points[5] = Point<2> (1,0.5);
  unit_points[6] = Point<2> (0.5,1);
  unit_points[7] = Point<2> (0,0.5);
  unit_points[8] = Point<2> (0.5,0.5);
};

  

template <>
void FEQuadraticSub<2>::get_support_points (const typename DoFHandler<2>::cell_iterator &cell,
					    const Boundary<2>&,
					    vector<Point<2> >  &support_points) const {
  Assert (support_points.size() == total_dofs,
	  ExcWrongFieldDimension (support_points.size(), total_dofs));
  
  for (unsigned int vertex=0; vertex<4; ++vertex)
    support_points[vertex] = cell->vertex(vertex);

				   // for the bilinear mapping, the centers
				   // of the face on the unit cell are mapped
				   // to the mean coordinates of the vertices
  for (unsigned int line=0; line<4; ++line)
    support_points[4+line] = (cell->line(line)->vertex(0) +
			      cell->line(line)->vertex(1)) / 2;
				   // same for the center of the square:
				   // since all four linear basis functions
				   // take on the value 1/4 at the center,
				   // the center is mapped to the mean
				   // coordinates of the four vertices
  support_points[8] = (support_points[0] +
		       support_points[1] +
		       support_points[2] +
		       support_points[3]) / 4;
};



template <>
void FEQuadraticSub<2>::get_face_support_points (const typename DoFHandler<2>::face_iterator &face,
						const Boundary<2>  &,
						vector<Point<2> >  &support_points) const {
  Assert (support_points.size() == dofs_per_face,
	  ExcWrongFieldDimension (support_points.size(), dofs_per_face));

  for (unsigned int vertex=0; vertex<2; ++vertex)
    support_points[vertex] = face->vertex(vertex);
  support_points[2] = (support_points[0] + support_points[1]) / 2;
};

#endif



#if deal_II_dimension == 3

template <>
FEQuadraticSub<3>::FEQuadraticSub () :
		FELinearMapping<3> (1, 1, 1, 1)
{
  interface_constraints(0,8) = 1.0;
  interface_constraints(1,4) = 1.0;
  interface_constraints(2,5) = 1.0;
  interface_constraints(3,6) = 1.0;
  interface_constraints(4,7) = 1.0;
  interface_constraints(5,4) = 3.0/8.0;
  interface_constraints(5,6) = -1.0/8.0;
  interface_constraints(5,8) = 3.0/4.0;
  interface_constraints(6,5) = 3.0/8.0;
  interface_constraints(6,7) = -1.0/8.0;
  interface_constraints(6,8) = 3.0/4.0;
  interface_constraints(7,4) = -1.0/8.0;
  interface_constraints(7,6) = 3.0/8.0;
  interface_constraints(7,8) = 3.0/4.0;
  interface_constraints(8,5) = -1.0/8.0;
  interface_constraints(8,7) = 3.0/8.0;
  interface_constraints(8,8) = 3.0/4.0;
  interface_constraints(9,0) = 3.0/8.0;
  interface_constraints(9,1) = -1.0/8.0;
  interface_constraints(9,4) = 3.0/4.0;
  interface_constraints(10,0) = -1.0/8.0;
  interface_constraints(10,1) = 3.0/8.0;
  interface_constraints(10,4) = 3.0/4.0;
  interface_constraints(11,1) = 3.0/8.0;
  interface_constraints(11,2) = -1.0/8.0;
  interface_constraints(11,5) = 3.0/4.0;
  interface_constraints(12,1) = -1.0/8.0;
  interface_constraints(12,2) = 3.0/8.0;
  interface_constraints(12,5) = 3.0/4.0;
  interface_constraints(13,2) = -1.0/8.0;
  interface_constraints(13,3) = 3.0/8.0;
  interface_constraints(13,6) = 3.0/4.0;
  interface_constraints(14,2) = 3.0/8.0;
  interface_constraints(14,3) = -1.0/8.0;
  interface_constraints(14,6) = 3.0/4.0;
  interface_constraints(15,0) = 3.0/8.0;
  interface_constraints(15,3) = -1.0/8.0;
  interface_constraints(15,7) = 3.0/4.0;
  interface_constraints(16,0) = -1.0/8.0;
  interface_constraints(16,3) = 3.0/8.0;
  interface_constraints(16,7) = 3.0/4.0;
  interface_constraints(17,0) = 9.0/64.0;
  interface_constraints(17,1) = -3.0/64.0;
  interface_constraints(17,2) = 1.0/64.0;
  interface_constraints(17,3) = -3.0/64.0;
  interface_constraints(17,4) = 9.0/32.0;
  interface_constraints(17,5) = -3.0/32.0;
  interface_constraints(17,6) = -3.0/32.0;
  interface_constraints(17,7) = 9.0/32.0;
  interface_constraints(17,8) = 9.0/16.0;
  interface_constraints(18,0) = -3.0/64.0;
  interface_constraints(18,1) = 9.0/64.0;
  interface_constraints(18,2) = -3.0/64.0;
  interface_constraints(18,3) = 1.0/64.0;
  interface_constraints(18,4) = 9.0/32.0;
  interface_constraints(18,5) = 9.0/32.0;
  interface_constraints(18,6) = -3.0/32.0;
  interface_constraints(18,7) = -3.0/32.0;
  interface_constraints(18,8) = 9.0/16.0;
  interface_constraints(19,0) = 1.0/64.0;
  interface_constraints(19,1) = -3.0/64.0;
  interface_constraints(19,2) = 9.0/64.0;
  interface_constraints(19,3) = -3.0/64.0;
  interface_constraints(19,4) = -3.0/32.0;
  interface_constraints(19,5) = 9.0/32.0;
  interface_constraints(19,6) = 9.0/32.0;
  interface_constraints(19,7) = -3.0/32.0;
  interface_constraints(19,8) = 9.0/16.0;
  interface_constraints(20,0) = -3.0/64.0;
  interface_constraints(20,1) = 1.0/64.0;
  interface_constraints(20,2) = -3.0/64.0;
  interface_constraints(20,3) = 9.0/64.0;
  interface_constraints(20,4) = -3.0/32.0;
  interface_constraints(20,5) = -3.0/32.0;
  interface_constraints(20,6) = 9.0/32.0;
  interface_constraints(20,7) = 9.0/32.0;
  interface_constraints(20,8) = 9.0/16.0;

  initialize_matrices ();
};



template <>
FEQuadraticSub<3>::FEQuadraticSub (const int) :
		FELinearMapping<3> (0, 0, 0, 27)
{
  initialize_matrices ();
};



template <>
void FEQuadraticSub<3>::initialize_matrices () {
  prolongation[0](0,0) = 1.0;
  prolongation[0](1,8) = 1.0;
  prolongation[0](2,20) = 1.0;
  prolongation[0](3,11) = 1.0;
  prolongation[0](4,16) = 1.0;
  prolongation[0](5,22) = 1.0;
  prolongation[0](6,26) = 1.0;
  prolongation[0](7,25) = 1.0;
  prolongation[0](8,0) = 3.0/8.0;
  prolongation[0](8,1) = -1.0/8.0;
  prolongation[0](8,8) = 3.0/4.0;
  prolongation[0](9,8) = 3.0/8.0;
  prolongation[0](9,10) = -1.0/8.0;
  prolongation[0](9,20) = 3.0/4.0;
  prolongation[0](10,9) = -1.0/8.0;
  prolongation[0](10,11) = 3.0/8.0;
  prolongation[0](10,20) = 3.0/4.0;
  prolongation[0](11,0) = 3.0/8.0;
  prolongation[0](11,3) = -1.0/8.0;
  prolongation[0](11,11) = 3.0/4.0;
  prolongation[0](12,16) = 3.0/8.0;
  prolongation[0](12,17) = -1.0/8.0;
  prolongation[0](12,22) = 3.0/4.0;
  prolongation[0](13,22) = 3.0/8.0;
  prolongation[0](13,24) = -1.0/8.0;
  prolongation[0](13,26) = 3.0/4.0;
  prolongation[0](14,23) = -1.0/8.0;
  prolongation[0](14,25) = 3.0/8.0;
  prolongation[0](14,26) = 3.0/4.0;
  prolongation[0](15,16) = 3.0/8.0;
  prolongation[0](15,19) = -1.0/8.0;
  prolongation[0](15,25) = 3.0/4.0;
  prolongation[0](16,0) = 3.0/8.0;
  prolongation[0](16,4) = -1.0/8.0;
  prolongation[0](16,16) = 3.0/4.0;
  prolongation[0](17,8) = 3.0/8.0;
  prolongation[0](17,12) = -1.0/8.0;
  prolongation[0](17,22) = 3.0/4.0;
  prolongation[0](18,20) = 3.0/8.0;
  prolongation[0](18,21) = -1.0/8.0;
  prolongation[0](18,26) = 3.0/4.0;
  prolongation[0](19,11) = 3.0/8.0;
  prolongation[0](19,15) = -1.0/8.0;
  prolongation[0](19,25) = 3.0/4.0;
  prolongation[0](20,0) = 9.0/64.0;
  prolongation[0](20,1) = -3.0/64.0;
  prolongation[0](20,2) = 1.0/64.0;
  prolongation[0](20,3) = -3.0/64.0;
  prolongation[0](20,8) = 9.0/32.0;
  prolongation[0](20,9) = -3.0/32.0;
  prolongation[0](20,10) = -3.0/32.0;
  prolongation[0](20,11) = 9.0/32.0;
  prolongation[0](20,20) = 9.0/16.0;
  prolongation[0](21,16) = 9.0/64.0;
  prolongation[0](21,17) = -3.0/64.0;
  prolongation[0](21,18) = 1.0/64.0;
  prolongation[0](21,19) = -3.0/64.0;
  prolongation[0](21,22) = 9.0/32.0;
  prolongation[0](21,23) = -3.0/32.0;
  prolongation[0](21,24) = -3.0/32.0;
  prolongation[0](21,25) = 9.0/32.0;
  prolongation[0](21,26) = 9.0/16.0;
  prolongation[0](22,0) = 9.0/64.0;
  prolongation[0](22,1) = -3.0/64.0;
  prolongation[0](22,4) = -3.0/64.0;
  prolongation[0](22,5) = 1.0/64.0;
  prolongation[0](22,8) = 9.0/32.0;
  prolongation[0](22,12) = -3.0/32.0;
  prolongation[0](22,16) = 9.0/32.0;
  prolongation[0](22,17) = -3.0/32.0;
  prolongation[0](22,22) = 9.0/16.0;
  prolongation[0](23,8) = 9.0/64.0;
  prolongation[0](23,10) = -3.0/64.0;
  prolongation[0](23,12) = -3.0/64.0;
  prolongation[0](23,14) = 1.0/64.0;
  prolongation[0](23,20) = 9.0/32.0;
  prolongation[0](23,21) = -3.0/32.0;
  prolongation[0](23,22) = 9.0/32.0;
  prolongation[0](23,24) = -3.0/32.0;
  prolongation[0](23,26) = 9.0/16.0;
  prolongation[0](24,9) = -3.0/64.0;
  prolongation[0](24,11) = 9.0/64.0;
  prolongation[0](24,13) = 1.0/64.0;
  prolongation[0](24,15) = -3.0/64.0;
  prolongation[0](24,20) = 9.0/32.0;
  prolongation[0](24,21) = -3.0/32.0;
  prolongation[0](24,23) = -3.0/32.0;
  prolongation[0](24,25) = 9.0/32.0;
  prolongation[0](24,26) = 9.0/16.0;
  prolongation[0](25,0) = 9.0/64.0;
  prolongation[0](25,3) = -3.0/64.0;
  prolongation[0](25,4) = -3.0/64.0;
  prolongation[0](25,7) = 1.0/64.0;
  prolongation[0](25,11) = 9.0/32.0;
  prolongation[0](25,15) = -3.0/32.0;
  prolongation[0](25,16) = 9.0/32.0;
  prolongation[0](25,19) = -3.0/32.0;
  prolongation[0](25,25) = 9.0/16.0;
  prolongation[0](26,0) = 27.0/512.0;
  prolongation[0](26,1) = -9.0/512.0;
  prolongation[0](26,2) = 3.0/512.0;
  prolongation[0](26,3) = -9.0/512.0;
  prolongation[0](26,4) = -9.0/512.0;
  prolongation[0](26,5) = 3.0/512.0;
  prolongation[0](26,6) = -1.0/512.0;
  prolongation[0](26,7) = 3.0/512.0;
  prolongation[0](26,8) = 27.0/256.0;
  prolongation[0](26,9) = -9.0/256.0;
  prolongation[0](26,10) = -9.0/256.0;
  prolongation[0](26,11) = 27.0/256.0;
  prolongation[0](26,12) = -9.0/256.0;
  prolongation[0](26,13) = 3.0/256.0;
  prolongation[0](26,14) = 3.0/256.0;
  prolongation[0](26,15) = -9.0/256.0;
  prolongation[0](26,16) = 27.0/256.0;
  prolongation[0](26,17) = -9.0/256.0;
  prolongation[0](26,18) = 3.0/256.0;
  prolongation[0](26,19) = -9.0/256.0;
  prolongation[0](26,20) = 27.0/128.0;
  prolongation[0](26,21) = -9.0/128.0;
  prolongation[0](26,22) = 27.0/128.0;
  prolongation[0](26,23) = -9.0/128.0;
  prolongation[0](26,24) = -9.0/128.0;
  prolongation[0](26,25) = 27.0/128.0;
  prolongation[0](26,26) = 27.0/64.0;
  prolongation[1](0,8) = 1.0;
  prolongation[1](1,1) = 1.0;
  prolongation[1](2,9) = 1.0;
  prolongation[1](3,20) = 1.0;
  prolongation[1](4,22) = 1.0;
  prolongation[1](5,17) = 1.0;
  prolongation[1](6,23) = 1.0;
  prolongation[1](7,26) = 1.0;
  prolongation[1](8,0) = -1.0/8.0;
  prolongation[1](8,1) = 3.0/8.0;
  prolongation[1](8,8) = 3.0/4.0;
  prolongation[1](9,1) = 3.0/8.0;
  prolongation[1](9,2) = -1.0/8.0;
  prolongation[1](9,9) = 3.0/4.0;
  prolongation[1](10,9) = 3.0/8.0;
  prolongation[1](10,11) = -1.0/8.0;
  prolongation[1](10,20) = 3.0/4.0;
  prolongation[1](11,8) = 3.0/8.0;
  prolongation[1](11,10) = -1.0/8.0;
  prolongation[1](11,20) = 3.0/4.0;
  prolongation[1](12,16) = -1.0/8.0;
  prolongation[1](12,17) = 3.0/8.0;
  prolongation[1](12,22) = 3.0/4.0;
  prolongation[1](13,17) = 3.0/8.0;
  prolongation[1](13,18) = -1.0/8.0;
  prolongation[1](13,23) = 3.0/4.0;
  prolongation[1](14,23) = 3.0/8.0;
  prolongation[1](14,25) = -1.0/8.0;
  prolongation[1](14,26) = 3.0/4.0;
  prolongation[1](15,22) = 3.0/8.0;
  prolongation[1](15,24) = -1.0/8.0;
  prolongation[1](15,26) = 3.0/4.0;
  prolongation[1](16,8) = 3.0/8.0;
  prolongation[1](16,12) = -1.0/8.0;
  prolongation[1](16,22) = 3.0/4.0;
  prolongation[1](17,1) = 3.0/8.0;
  prolongation[1](17,5) = -1.0/8.0;
  prolongation[1](17,17) = 3.0/4.0;
  prolongation[1](18,9) = 3.0/8.0;
  prolongation[1](18,13) = -1.0/8.0;
  prolongation[1](18,23) = 3.0/4.0;
  prolongation[1](19,20) = 3.0/8.0;
  prolongation[1](19,21) = -1.0/8.0;
  prolongation[1](19,26) = 3.0/4.0;
  prolongation[1](20,0) = -3.0/64.0;
  prolongation[1](20,1) = 9.0/64.0;
  prolongation[1](20,2) = -3.0/64.0;
  prolongation[1](20,3) = 1.0/64.0;
  prolongation[1](20,8) = 9.0/32.0;
  prolongation[1](20,9) = 9.0/32.0;
  prolongation[1](20,10) = -3.0/32.0;
  prolongation[1](20,11) = -3.0/32.0;
  prolongation[1](20,20) = 9.0/16.0;
  prolongation[1](21,16) = -3.0/64.0;
  prolongation[1](21,17) = 9.0/64.0;
  prolongation[1](21,18) = -3.0/64.0;
  prolongation[1](21,19) = 1.0/64.0;
  prolongation[1](21,22) = 9.0/32.0;
  prolongation[1](21,23) = 9.0/32.0;
  prolongation[1](21,24) = -3.0/32.0;
  prolongation[1](21,25) = -3.0/32.0;
  prolongation[1](21,26) = 9.0/16.0;
  prolongation[1](22,0) = -3.0/64.0;
  prolongation[1](22,1) = 9.0/64.0;
  prolongation[1](22,4) = 1.0/64.0;
  prolongation[1](22,5) = -3.0/64.0;
  prolongation[1](22,8) = 9.0/32.0;
  prolongation[1](22,12) = -3.0/32.0;
  prolongation[1](22,16) = -3.0/32.0;
  prolongation[1](22,17) = 9.0/32.0;
  prolongation[1](22,22) = 9.0/16.0;
  prolongation[1](23,1) = 9.0/64.0;
  prolongation[1](23,2) = -3.0/64.0;
  prolongation[1](23,5) = -3.0/64.0;
  prolongation[1](23,6) = 1.0/64.0;
  prolongation[1](23,9) = 9.0/32.0;
  prolongation[1](23,13) = -3.0/32.0;
  prolongation[1](23,17) = 9.0/32.0;
  prolongation[1](23,18) = -3.0/32.0;
  prolongation[1](23,23) = 9.0/16.0;
  prolongation[1](24,9) = 9.0/64.0;
  prolongation[1](24,11) = -3.0/64.0;
  prolongation[1](24,13) = -3.0/64.0;
  prolongation[1](24,15) = 1.0/64.0;
  prolongation[1](24,20) = 9.0/32.0;
  prolongation[1](24,21) = -3.0/32.0;
  prolongation[1](24,23) = 9.0/32.0;
  prolongation[1](24,25) = -3.0/32.0;
  prolongation[1](24,26) = 9.0/16.0;
  prolongation[1](25,8) = 9.0/64.0;
  prolongation[1](25,10) = -3.0/64.0;
  prolongation[1](25,12) = -3.0/64.0;
  prolongation[1](25,14) = 1.0/64.0;
  prolongation[1](25,20) = 9.0/32.0;
  prolongation[1](25,21) = -3.0/32.0;
  prolongation[1](25,22) = 9.0/32.0;
  prolongation[1](25,24) = -3.0/32.0;
  prolongation[1](25,26) = 9.0/16.0;
  prolongation[1](26,0) = -9.0/512.0;
  prolongation[1](26,1) = 27.0/512.0;
  prolongation[1](26,2) = -9.0/512.0;
  prolongation[1](26,3) = 3.0/512.0;
  prolongation[1](26,4) = 3.0/512.0;
  prolongation[1](26,5) = -9.0/512.0;
  prolongation[1](26,6) = 3.0/512.0;
  prolongation[1](26,7) = -1.0/512.0;
  prolongation[1](26,8) = 27.0/256.0;
  prolongation[1](26,9) = 27.0/256.0;
  prolongation[1](26,10) = -9.0/256.0;
  prolongation[1](26,11) = -9.0/256.0;
  prolongation[1](26,12) = -9.0/256.0;
  prolongation[1](26,13) = -9.0/256.0;
  prolongation[1](26,14) = 3.0/256.0;
  prolongation[1](26,15) = 3.0/256.0;
  prolongation[1](26,16) = -9.0/256.0;
  prolongation[1](26,17) = 27.0/256.0;
  prolongation[1](26,18) = -9.0/256.0;
  prolongation[1](26,19) = 3.0/256.0;
  prolongation[1](26,20) = 27.0/128.0;
  prolongation[1](26,21) = -9.0/128.0;
  prolongation[1](26,22) = 27.0/128.0;
  prolongation[1](26,23) = 27.0/128.0;
  prolongation[1](26,24) = -9.0/128.0;
  prolongation[1](26,25) = -9.0/128.0;
  prolongation[1](26,26) = 27.0/64.0;
  prolongation[2](0,20) = 1.0;
  prolongation[2](1,9) = 1.0;
  prolongation[2](2,2) = 1.0;
  prolongation[2](3,10) = 1.0;
  prolongation[2](4,26) = 1.0;
  prolongation[2](5,23) = 1.0;
  prolongation[2](6,18) = 1.0;
  prolongation[2](7,24) = 1.0;
  prolongation[2](8,9) = 3.0/8.0;
  prolongation[2](8,11) = -1.0/8.0;
  prolongation[2](8,20) = 3.0/4.0;
  prolongation[2](9,1) = -1.0/8.0;
  prolongation[2](9,2) = 3.0/8.0;
  prolongation[2](9,9) = 3.0/4.0;
  prolongation[2](10,2) = 3.0/8.0;
  prolongation[2](10,3) = -1.0/8.0;
  prolongation[2](10,10) = 3.0/4.0;
  prolongation[2](11,8) = -1.0/8.0;
  prolongation[2](11,10) = 3.0/8.0;
  prolongation[2](11,20) = 3.0/4.0;
  prolongation[2](12,23) = 3.0/8.0;
  prolongation[2](12,25) = -1.0/8.0;
  prolongation[2](12,26) = 3.0/4.0;
  prolongation[2](13,17) = -1.0/8.0;
  prolongation[2](13,18) = 3.0/8.0;
  prolongation[2](13,23) = 3.0/4.0;
  prolongation[2](14,18) = 3.0/8.0;
  prolongation[2](14,19) = -1.0/8.0;
  prolongation[2](14,24) = 3.0/4.0;
  prolongation[2](15,22) = -1.0/8.0;
  prolongation[2](15,24) = 3.0/8.0;
  prolongation[2](15,26) = 3.0/4.0;
  prolongation[2](16,20) = 3.0/8.0;
  prolongation[2](16,21) = -1.0/8.0;
  prolongation[2](16,26) = 3.0/4.0;
  prolongation[2](17,9) = 3.0/8.0;
  prolongation[2](17,13) = -1.0/8.0;
  prolongation[2](17,23) = 3.0/4.0;
  prolongation[2](18,2) = 3.0/8.0;
  prolongation[2](18,6) = -1.0/8.0;
  prolongation[2](18,18) = 3.0/4.0;
  prolongation[2](19,10) = 3.0/8.0;
  prolongation[2](19,14) = -1.0/8.0;
  prolongation[2](19,24) = 3.0/4.0;
  prolongation[2](20,0) = 1.0/64.0;
  prolongation[2](20,1) = -3.0/64.0;
  prolongation[2](20,2) = 9.0/64.0;
  prolongation[2](20,3) = -3.0/64.0;
  prolongation[2](20,8) = -3.0/32.0;
  prolongation[2](20,9) = 9.0/32.0;
  prolongation[2](20,10) = 9.0/32.0;
  prolongation[2](20,11) = -3.0/32.0;
  prolongation[2](20,20) = 9.0/16.0;
  prolongation[2](21,16) = 1.0/64.0;
  prolongation[2](21,17) = -3.0/64.0;
  prolongation[2](21,18) = 9.0/64.0;
  prolongation[2](21,19) = -3.0/64.0;
  prolongation[2](21,22) = -3.0/32.0;
  prolongation[2](21,23) = 9.0/32.0;
  prolongation[2](21,24) = 9.0/32.0;
  prolongation[2](21,25) = -3.0/32.0;
  prolongation[2](21,26) = 9.0/16.0;
  prolongation[2](22,9) = 9.0/64.0;
  prolongation[2](22,11) = -3.0/64.0;
  prolongation[2](22,13) = -3.0/64.0;
  prolongation[2](22,15) = 1.0/64.0;
  prolongation[2](22,20) = 9.0/32.0;
  prolongation[2](22,21) = -3.0/32.0;
  prolongation[2](22,23) = 9.0/32.0;
  prolongation[2](22,25) = -3.0/32.0;
  prolongation[2](22,26) = 9.0/16.0;
  prolongation[2](23,1) = -3.0/64.0;
  prolongation[2](23,2) = 9.0/64.0;
  prolongation[2](23,5) = 1.0/64.0;
  prolongation[2](23,6) = -3.0/64.0;
  prolongation[2](23,9) = 9.0/32.0;
  prolongation[2](23,13) = -3.0/32.0;
  prolongation[2](23,17) = -3.0/32.0;
  prolongation[2](23,18) = 9.0/32.0;
  prolongation[2](23,23) = 9.0/16.0;
  prolongation[2](24,2) = 9.0/64.0;
  prolongation[2](24,3) = -3.0/64.0;
  prolongation[2](24,6) = -3.0/64.0;
  prolongation[2](24,7) = 1.0/64.0;
  prolongation[2](24,10) = 9.0/32.0;
  prolongation[2](24,14) = -3.0/32.0;
  prolongation[2](24,18) = 9.0/32.0;
  prolongation[2](24,19) = -3.0/32.0;
  prolongation[2](24,24) = 9.0/16.0;
  prolongation[2](25,8) = -3.0/64.0;
  prolongation[2](25,10) = 9.0/64.0;
  prolongation[2](25,12) = 1.0/64.0;
  prolongation[2](25,14) = -3.0/64.0;
  prolongation[2](25,20) = 9.0/32.0;
  prolongation[2](25,21) = -3.0/32.0;
  prolongation[2](25,22) = -3.0/32.0;
  prolongation[2](25,24) = 9.0/32.0;
  prolongation[2](25,26) = 9.0/16.0;
  prolongation[2](26,0) = 3.0/512.0;
  prolongation[2](26,1) = -9.0/512.0;
  prolongation[2](26,2) = 27.0/512.0;
  prolongation[2](26,3) = -9.0/512.0;
  prolongation[2](26,4) = -1.0/512.0;
  prolongation[2](26,5) = 3.0/512.0;
  prolongation[2](26,6) = -9.0/512.0;
  prolongation[2](26,7) = 3.0/512.0;
  prolongation[2](26,8) = -9.0/256.0;
  prolongation[2](26,9) = 27.0/256.0;
  prolongation[2](26,10) = 27.0/256.0;
  prolongation[2](26,11) = -9.0/256.0;
  prolongation[2](26,12) = 3.0/256.0;
  prolongation[2](26,13) = -9.0/256.0;
  prolongation[2](26,14) = -9.0/256.0;
  prolongation[2](26,15) = 3.0/256.0;
  prolongation[2](26,16) = 3.0/256.0;
  prolongation[2](26,17) = -9.0/256.0;
  prolongation[2](26,18) = 27.0/256.0;
  prolongation[2](26,19) = -9.0/256.0;
  prolongation[2](26,20) = 27.0/128.0;
  prolongation[2](26,21) = -9.0/128.0;
  prolongation[2](26,22) = -9.0/128.0;
  prolongation[2](26,23) = 27.0/128.0;
  prolongation[2](26,24) = 27.0/128.0;
  prolongation[2](26,25) = -9.0/128.0;
  prolongation[2](26,26) = 27.0/64.0;
  prolongation[3](0,11) = 1.0;
  prolongation[3](1,20) = 1.0;
  prolongation[3](2,10) = 1.0;
  prolongation[3](3,3) = 1.0;
  prolongation[3](4,25) = 1.0;
  prolongation[3](5,26) = 1.0;
  prolongation[3](6,24) = 1.0;
  prolongation[3](7,19) = 1.0;
  prolongation[3](8,9) = -1.0/8.0;
  prolongation[3](8,11) = 3.0/8.0;
  prolongation[3](8,20) = 3.0/4.0;
  prolongation[3](9,8) = -1.0/8.0;
  prolongation[3](9,10) = 3.0/8.0;
  prolongation[3](9,20) = 3.0/4.0;
  prolongation[3](10,2) = -1.0/8.0;
  prolongation[3](10,3) = 3.0/8.0;
  prolongation[3](10,10) = 3.0/4.0;
  prolongation[3](11,0) = -1.0/8.0;
  prolongation[3](11,3) = 3.0/8.0;
  prolongation[3](11,11) = 3.0/4.0;
  prolongation[3](12,23) = -1.0/8.0;
  prolongation[3](12,25) = 3.0/8.0;
  prolongation[3](12,26) = 3.0/4.0;
  prolongation[3](13,22) = -1.0/8.0;
  prolongation[3](13,24) = 3.0/8.0;
  prolongation[3](13,26) = 3.0/4.0;
  prolongation[3](14,18) = -1.0/8.0;
  prolongation[3](14,19) = 3.0/8.0;
  prolongation[3](14,24) = 3.0/4.0;
  prolongation[3](15,16) = -1.0/8.0;
  prolongation[3](15,19) = 3.0/8.0;
  prolongation[3](15,25) = 3.0/4.0;
  prolongation[3](16,11) = 3.0/8.0;
  prolongation[3](16,15) = -1.0/8.0;
  prolongation[3](16,25) = 3.0/4.0;
  prolongation[3](17,20) = 3.0/8.0;
  prolongation[3](17,21) = -1.0/8.0;
  prolongation[3](17,26) = 3.0/4.0;
  prolongation[3](18,10) = 3.0/8.0;
  prolongation[3](18,14) = -1.0/8.0;
  prolongation[3](18,24) = 3.0/4.0;
  prolongation[3](19,3) = 3.0/8.0;
  prolongation[3](19,7) = -1.0/8.0;
  prolongation[3](19,19) = 3.0/4.0;
  prolongation[3](20,0) = -3.0/64.0;
  prolongation[3](20,1) = 1.0/64.0;
  prolongation[3](20,2) = -3.0/64.0;
  prolongation[3](20,3) = 9.0/64.0;
  prolongation[3](20,8) = -3.0/32.0;
  prolongation[3](20,9) = -3.0/32.0;
  prolongation[3](20,10) = 9.0/32.0;
  prolongation[3](20,11) = 9.0/32.0;
  prolongation[3](20,20) = 9.0/16.0;
  prolongation[3](21,16) = -3.0/64.0;
  prolongation[3](21,17) = 1.0/64.0;
  prolongation[3](21,18) = -3.0/64.0;
  prolongation[3](21,19) = 9.0/64.0;
  prolongation[3](21,22) = -3.0/32.0;
  prolongation[3](21,23) = -3.0/32.0;
  prolongation[3](21,24) = 9.0/32.0;
  prolongation[3](21,25) = 9.0/32.0;
  prolongation[3](21,26) = 9.0/16.0;
  prolongation[3](22,9) = -3.0/64.0;
  prolongation[3](22,11) = 9.0/64.0;
  prolongation[3](22,13) = 1.0/64.0;
  prolongation[3](22,15) = -3.0/64.0;
  prolongation[3](22,20) = 9.0/32.0;
  prolongation[3](22,21) = -3.0/32.0;
  prolongation[3](22,23) = -3.0/32.0;
  prolongation[3](22,25) = 9.0/32.0;
  prolongation[3](22,26) = 9.0/16.0;
  prolongation[3](23,8) = -3.0/64.0;
  prolongation[3](23,10) = 9.0/64.0;
  prolongation[3](23,12) = 1.0/64.0;
  prolongation[3](23,14) = -3.0/64.0;
  prolongation[3](23,20) = 9.0/32.0;
  prolongation[3](23,21) = -3.0/32.0;
  prolongation[3](23,22) = -3.0/32.0;
  prolongation[3](23,24) = 9.0/32.0;
  prolongation[3](23,26) = 9.0/16.0;
  prolongation[3](24,2) = -3.0/64.0;
  prolongation[3](24,3) = 9.0/64.0;
  prolongation[3](24,6) = 1.0/64.0;
  prolongation[3](24,7) = -3.0/64.0;
  prolongation[3](24,10) = 9.0/32.0;
  prolongation[3](24,14) = -3.0/32.0;
  prolongation[3](24,18) = -3.0/32.0;
  prolongation[3](24,19) = 9.0/32.0;
  prolongation[3](24,24) = 9.0/16.0;
  prolongation[3](25,0) = -3.0/64.0;
  prolongation[3](25,3) = 9.0/64.0;
  prolongation[3](25,4) = 1.0/64.0;
  prolongation[3](25,7) = -3.0/64.0;
  prolongation[3](25,11) = 9.0/32.0;
  prolongation[3](25,15) = -3.0/32.0;
  prolongation[3](25,16) = -3.0/32.0;
  prolongation[3](25,19) = 9.0/32.0;
  prolongation[3](25,25) = 9.0/16.0;
  prolongation[3](26,0) = -9.0/512.0;
  prolongation[3](26,1) = 3.0/512.0;
  prolongation[3](26,2) = -9.0/512.0;
  prolongation[3](26,3) = 27.0/512.0;
  prolongation[3](26,4) = 3.0/512.0;
  prolongation[3](26,5) = -1.0/512.0;
  prolongation[3](26,6) = 3.0/512.0;
  prolongation[3](26,7) = -9.0/512.0;
  prolongation[3](26,8) = -9.0/256.0;
  prolongation[3](26,9) = -9.0/256.0;
  prolongation[3](26,10) = 27.0/256.0;
  prolongation[3](26,11) = 27.0/256.0;
  prolongation[3](26,12) = 3.0/256.0;
  prolongation[3](26,13) = 3.0/256.0;
  prolongation[3](26,14) = -9.0/256.0;
  prolongation[3](26,15) = -9.0/256.0;
  prolongation[3](26,16) = -9.0/256.0;
  prolongation[3](26,17) = 3.0/256.0;
  prolongation[3](26,18) = -9.0/256.0;
  prolongation[3](26,19) = 27.0/256.0;
  prolongation[3](26,20) = 27.0/128.0;
  prolongation[3](26,21) = -9.0/128.0;
  prolongation[3](26,22) = -9.0/128.0;
  prolongation[3](26,23) = -9.0/128.0;
  prolongation[3](26,24) = 27.0/128.0;
  prolongation[3](26,25) = 27.0/128.0;
  prolongation[3](26,26) = 27.0/64.0;
  prolongation[4](0,16) = 1.0;
  prolongation[4](1,22) = 1.0;
  prolongation[4](2,26) = 1.0;
  prolongation[4](3,25) = 1.0;
  prolongation[4](4,4) = 1.0;
  prolongation[4](5,12) = 1.0;
  prolongation[4](6,21) = 1.0;
  prolongation[4](7,15) = 1.0;
  prolongation[4](8,16) = 3.0/8.0;
  prolongation[4](8,17) = -1.0/8.0;
  prolongation[4](8,22) = 3.0/4.0;
  prolongation[4](9,22) = 3.0/8.0;
  prolongation[4](9,24) = -1.0/8.0;
  prolongation[4](9,26) = 3.0/4.0;
  prolongation[4](10,23) = -1.0/8.0;
  prolongation[4](10,25) = 3.0/8.0;
  prolongation[4](10,26) = 3.0/4.0;
  prolongation[4](11,16) = 3.0/8.0;
  prolongation[4](11,19) = -1.0/8.0;
  prolongation[4](11,25) = 3.0/4.0;
  prolongation[4](12,4) = 3.0/8.0;
  prolongation[4](12,5) = -1.0/8.0;
  prolongation[4](12,12) = 3.0/4.0;
  prolongation[4](13,12) = 3.0/8.0;
  prolongation[4](13,14) = -1.0/8.0;
  prolongation[4](13,21) = 3.0/4.0;
  prolongation[4](14,13) = -1.0/8.0;
  prolongation[4](14,15) = 3.0/8.0;
  prolongation[4](14,21) = 3.0/4.0;
  prolongation[4](15,4) = 3.0/8.0;
  prolongation[4](15,7) = -1.0/8.0;
  prolongation[4](15,15) = 3.0/4.0;
  prolongation[4](16,0) = -1.0/8.0;
  prolongation[4](16,4) = 3.0/8.0;
  prolongation[4](16,16) = 3.0/4.0;
  prolongation[4](17,8) = -1.0/8.0;
  prolongation[4](17,12) = 3.0/8.0;
  prolongation[4](17,22) = 3.0/4.0;
  prolongation[4](18,20) = -1.0/8.0;
  prolongation[4](18,21) = 3.0/8.0;
  prolongation[4](18,26) = 3.0/4.0;
  prolongation[4](19,11) = -1.0/8.0;
  prolongation[4](19,15) = 3.0/8.0;
  prolongation[4](19,25) = 3.0/4.0;
  prolongation[4](20,16) = 9.0/64.0;
  prolongation[4](20,17) = -3.0/64.0;
  prolongation[4](20,18) = 1.0/64.0;
  prolongation[4](20,19) = -3.0/64.0;
  prolongation[4](20,22) = 9.0/32.0;
  prolongation[4](20,23) = -3.0/32.0;
  prolongation[4](20,24) = -3.0/32.0;
  prolongation[4](20,25) = 9.0/32.0;
  prolongation[4](20,26) = 9.0/16.0;
  prolongation[4](21,4) = 9.0/64.0;
  prolongation[4](21,5) = -3.0/64.0;
  prolongation[4](21,6) = 1.0/64.0;
  prolongation[4](21,7) = -3.0/64.0;
  prolongation[4](21,12) = 9.0/32.0;
  prolongation[4](21,13) = -3.0/32.0;
  prolongation[4](21,14) = -3.0/32.0;
  prolongation[4](21,15) = 9.0/32.0;
  prolongation[4](21,21) = 9.0/16.0;
  prolongation[4](22,0) = -3.0/64.0;
  prolongation[4](22,1) = 1.0/64.0;
  prolongation[4](22,4) = 9.0/64.0;
  prolongation[4](22,5) = -3.0/64.0;
  prolongation[4](22,8) = -3.0/32.0;
  prolongation[4](22,12) = 9.0/32.0;
  prolongation[4](22,16) = 9.0/32.0;
  prolongation[4](22,17) = -3.0/32.0;
  prolongation[4](22,22) = 9.0/16.0;
  prolongation[4](23,8) = -3.0/64.0;
  prolongation[4](23,10) = 1.0/64.0;
  prolongation[4](23,12) = 9.0/64.0;
  prolongation[4](23,14) = -3.0/64.0;
  prolongation[4](23,20) = -3.0/32.0;
  prolongation[4](23,21) = 9.0/32.0;
  prolongation[4](23,22) = 9.0/32.0;
  prolongation[4](23,24) = -3.0/32.0;
  prolongation[4](23,26) = 9.0/16.0;
  prolongation[4](24,9) = 1.0/64.0;
  prolongation[4](24,11) = -3.0/64.0;
  prolongation[4](24,13) = -3.0/64.0;
  prolongation[4](24,15) = 9.0/64.0;
  prolongation[4](24,20) = -3.0/32.0;
  prolongation[4](24,21) = 9.0/32.0;
  prolongation[4](24,23) = -3.0/32.0;
  prolongation[4](24,25) = 9.0/32.0;
  prolongation[4](24,26) = 9.0/16.0;
  prolongation[4](25,0) = -3.0/64.0;
  prolongation[4](25,3) = 1.0/64.0;
  prolongation[4](25,4) = 9.0/64.0;
  prolongation[4](25,7) = -3.0/64.0;
  prolongation[4](25,11) = -3.0/32.0;
  prolongation[4](25,15) = 9.0/32.0;
  prolongation[4](25,16) = 9.0/32.0;
  prolongation[4](25,19) = -3.0/32.0;
  prolongation[4](25,25) = 9.0/16.0;
  prolongation[4](26,0) = -9.0/512.0;
  prolongation[4](26,1) = 3.0/512.0;
  prolongation[4](26,2) = -1.0/512.0;
  prolongation[4](26,3) = 3.0/512.0;
  prolongation[4](26,4) = 27.0/512.0;
  prolongation[4](26,5) = -9.0/512.0;
  prolongation[4](26,6) = 3.0/512.0;
  prolongation[4](26,7) = -9.0/512.0;
  prolongation[4](26,8) = -9.0/256.0;
  prolongation[4](26,9) = 3.0/256.0;
  prolongation[4](26,10) = 3.0/256.0;
  prolongation[4](26,11) = -9.0/256.0;
  prolongation[4](26,12) = 27.0/256.0;
  prolongation[4](26,13) = -9.0/256.0;
  prolongation[4](26,14) = -9.0/256.0;
  prolongation[4](26,15) = 27.0/256.0;
  prolongation[4](26,16) = 27.0/256.0;
  prolongation[4](26,17) = -9.0/256.0;
  prolongation[4](26,18) = 3.0/256.0;
  prolongation[4](26,19) = -9.0/256.0;
  prolongation[4](26,20) = -9.0/128.0;
  prolongation[4](26,21) = 27.0/128.0;
  prolongation[4](26,22) = 27.0/128.0;
  prolongation[4](26,23) = -9.0/128.0;
  prolongation[4](26,24) = -9.0/128.0;
  prolongation[4](26,25) = 27.0/128.0;
  prolongation[4](26,26) = 27.0/64.0;
  prolongation[5](0,22) = 1.0;
  prolongation[5](1,17) = 1.0;
  prolongation[5](2,23) = 1.0;
  prolongation[5](3,26) = 1.0;
  prolongation[5](4,12) = 1.0;
  prolongation[5](5,5) = 1.0;
  prolongation[5](6,13) = 1.0;
  prolongation[5](7,21) = 1.0;
  prolongation[5](8,16) = -1.0/8.0;
  prolongation[5](8,17) = 3.0/8.0;
  prolongation[5](8,22) = 3.0/4.0;
  prolongation[5](9,17) = 3.0/8.0;
  prolongation[5](9,18) = -1.0/8.0;
  prolongation[5](9,23) = 3.0/4.0;
  prolongation[5](10,23) = 3.0/8.0;
  prolongation[5](10,25) = -1.0/8.0;
  prolongation[5](10,26) = 3.0/4.0;
  prolongation[5](11,22) = 3.0/8.0;
  prolongation[5](11,24) = -1.0/8.0;
  prolongation[5](11,26) = 3.0/4.0;
  prolongation[5](12,4) = -1.0/8.0;
  prolongation[5](12,5) = 3.0/8.0;
  prolongation[5](12,12) = 3.0/4.0;
  prolongation[5](13,5) = 3.0/8.0;
  prolongation[5](13,6) = -1.0/8.0;
  prolongation[5](13,13) = 3.0/4.0;
  prolongation[5](14,13) = 3.0/8.0;
  prolongation[5](14,15) = -1.0/8.0;
  prolongation[5](14,21) = 3.0/4.0;
  prolongation[5](15,12) = 3.0/8.0;
  prolongation[5](15,14) = -1.0/8.0;
  prolongation[5](15,21) = 3.0/4.0;
  prolongation[5](16,8) = -1.0/8.0;
  prolongation[5](16,12) = 3.0/8.0;
  prolongation[5](16,22) = 3.0/4.0;
  prolongation[5](17,1) = -1.0/8.0;
  prolongation[5](17,5) = 3.0/8.0;
  prolongation[5](17,17) = 3.0/4.0;
  prolongation[5](18,9) = -1.0/8.0;
  prolongation[5](18,13) = 3.0/8.0;
  prolongation[5](18,23) = 3.0/4.0;
  prolongation[5](19,20) = -1.0/8.0;
  prolongation[5](19,21) = 3.0/8.0;
  prolongation[5](19,26) = 3.0/4.0;
  prolongation[5](20,16) = -3.0/64.0;
  prolongation[5](20,17) = 9.0/64.0;
  prolongation[5](20,18) = -3.0/64.0;
  prolongation[5](20,19) = 1.0/64.0;
  prolongation[5](20,22) = 9.0/32.0;
  prolongation[5](20,23) = 9.0/32.0;
  prolongation[5](20,24) = -3.0/32.0;
  prolongation[5](20,25) = -3.0/32.0;
  prolongation[5](20,26) = 9.0/16.0;
  prolongation[5](21,4) = -3.0/64.0;
  prolongation[5](21,5) = 9.0/64.0;
  prolongation[5](21,6) = -3.0/64.0;
  prolongation[5](21,7) = 1.0/64.0;
  prolongation[5](21,12) = 9.0/32.0;
  prolongation[5](21,13) = 9.0/32.0;
  prolongation[5](21,14) = -3.0/32.0;
  prolongation[5](21,15) = -3.0/32.0;
  prolongation[5](21,21) = 9.0/16.0;
  prolongation[5](22,0) = 1.0/64.0;
  prolongation[5](22,1) = -3.0/64.0;
  prolongation[5](22,4) = -3.0/64.0;
  prolongation[5](22,5) = 9.0/64.0;
  prolongation[5](22,8) = -3.0/32.0;
  prolongation[5](22,12) = 9.0/32.0;
  prolongation[5](22,16) = -3.0/32.0;
  prolongation[5](22,17) = 9.0/32.0;
  prolongation[5](22,22) = 9.0/16.0;
  prolongation[5](23,1) = -3.0/64.0;
  prolongation[5](23,2) = 1.0/64.0;
  prolongation[5](23,5) = 9.0/64.0;
  prolongation[5](23,6) = -3.0/64.0;
  prolongation[5](23,9) = -3.0/32.0;
  prolongation[5](23,13) = 9.0/32.0;
  prolongation[5](23,17) = 9.0/32.0;
  prolongation[5](23,18) = -3.0/32.0;
  prolongation[5](23,23) = 9.0/16.0;
  prolongation[5](24,9) = -3.0/64.0;
  prolongation[5](24,11) = 1.0/64.0;
  prolongation[5](24,13) = 9.0/64.0;
  prolongation[5](24,15) = -3.0/64.0;
  prolongation[5](24,20) = -3.0/32.0;
  prolongation[5](24,21) = 9.0/32.0;
  prolongation[5](24,23) = 9.0/32.0;
  prolongation[5](24,25) = -3.0/32.0;
  prolongation[5](24,26) = 9.0/16.0;
  prolongation[5](25,8) = -3.0/64.0;
  prolongation[5](25,10) = 1.0/64.0;
  prolongation[5](25,12) = 9.0/64.0;
  prolongation[5](25,14) = -3.0/64.0;
  prolongation[5](25,20) = -3.0/32.0;
  prolongation[5](25,21) = 9.0/32.0;
  prolongation[5](25,22) = 9.0/32.0;
  prolongation[5](25,24) = -3.0/32.0;
  prolongation[5](25,26) = 9.0/16.0;
  prolongation[5](26,0) = 3.0/512.0;
  prolongation[5](26,1) = -9.0/512.0;
  prolongation[5](26,2) = 3.0/512.0;
  prolongation[5](26,3) = -1.0/512.0;
  prolongation[5](26,4) = -9.0/512.0;
  prolongation[5](26,5) = 27.0/512.0;
  prolongation[5](26,6) = -9.0/512.0;
  prolongation[5](26,7) = 3.0/512.0;
  prolongation[5](26,8) = -9.0/256.0;
  prolongation[5](26,9) = -9.0/256.0;
  prolongation[5](26,10) = 3.0/256.0;
  prolongation[5](26,11) = 3.0/256.0;
  prolongation[5](26,12) = 27.0/256.0;
  prolongation[5](26,13) = 27.0/256.0;
  prolongation[5](26,14) = -9.0/256.0;
  prolongation[5](26,15) = -9.0/256.0;
  prolongation[5](26,16) = -9.0/256.0;
  prolongation[5](26,17) = 27.0/256.0;
  prolongation[5](26,18) = -9.0/256.0;
  prolongation[5](26,19) = 3.0/256.0;
  prolongation[5](26,20) = -9.0/128.0;
  prolongation[5](26,21) = 27.0/128.0;
  prolongation[5](26,22) = 27.0/128.0;
  prolongation[5](26,23) = 27.0/128.0;
  prolongation[5](26,24) = -9.0/128.0;
  prolongation[5](26,25) = -9.0/128.0;
  prolongation[5](26,26) = 27.0/64.0;
  prolongation[6](0,26) = 1.0;
  prolongation[6](1,23) = 1.0;
  prolongation[6](2,18) = 1.0;
  prolongation[6](3,24) = 1.0;
  prolongation[6](4,21) = 1.0;
  prolongation[6](5,13) = 1.0;
  prolongation[6](6,6) = 1.0;
  prolongation[6](7,14) = 1.0;
  prolongation[6](8,23) = 3.0/8.0;
  prolongation[6](8,25) = -1.0/8.0;
  prolongation[6](8,26) = 3.0/4.0;
  prolongation[6](9,17) = -1.0/8.0;
  prolongation[6](9,18) = 3.0/8.0;
  prolongation[6](9,23) = 3.0/4.0;
  prolongation[6](10,18) = 3.0/8.0;
  prolongation[6](10,19) = -1.0/8.0;
  prolongation[6](10,24) = 3.0/4.0;
  prolongation[6](11,22) = -1.0/8.0;
  prolongation[6](11,24) = 3.0/8.0;
  prolongation[6](11,26) = 3.0/4.0;
  prolongation[6](12,13) = 3.0/8.0;
  prolongation[6](12,15) = -1.0/8.0;
  prolongation[6](12,21) = 3.0/4.0;
  prolongation[6](13,5) = -1.0/8.0;
  prolongation[6](13,6) = 3.0/8.0;
  prolongation[6](13,13) = 3.0/4.0;
  prolongation[6](14,6) = 3.0/8.0;
  prolongation[6](14,7) = -1.0/8.0;
  prolongation[6](14,14) = 3.0/4.0;
  prolongation[6](15,12) = -1.0/8.0;
  prolongation[6](15,14) = 3.0/8.0;
  prolongation[6](15,21) = 3.0/4.0;
  prolongation[6](16,20) = -1.0/8.0;
  prolongation[6](16,21) = 3.0/8.0;
  prolongation[6](16,26) = 3.0/4.0;
  prolongation[6](17,9) = -1.0/8.0;
  prolongation[6](17,13) = 3.0/8.0;
  prolongation[6](17,23) = 3.0/4.0;
  prolongation[6](18,2) = -1.0/8.0;
  prolongation[6](18,6) = 3.0/8.0;
  prolongation[6](18,18) = 3.0/4.0;
  prolongation[6](19,10) = -1.0/8.0;
  prolongation[6](19,14) = 3.0/8.0;
  prolongation[6](19,24) = 3.0/4.0;
  prolongation[6](20,16) = 1.0/64.0;
  prolongation[6](20,17) = -3.0/64.0;
  prolongation[6](20,18) = 9.0/64.0;
  prolongation[6](20,19) = -3.0/64.0;
  prolongation[6](20,22) = -3.0/32.0;
  prolongation[6](20,23) = 9.0/32.0;
  prolongation[6](20,24) = 9.0/32.0;
  prolongation[6](20,25) = -3.0/32.0;
  prolongation[6](20,26) = 9.0/16.0;
  prolongation[6](21,4) = 1.0/64.0;
  prolongation[6](21,5) = -3.0/64.0;
  prolongation[6](21,6) = 9.0/64.0;
  prolongation[6](21,7) = -3.0/64.0;
  prolongation[6](21,12) = -3.0/32.0;
  prolongation[6](21,13) = 9.0/32.0;
  prolongation[6](21,14) = 9.0/32.0;
  prolongation[6](21,15) = -3.0/32.0;
  prolongation[6](21,21) = 9.0/16.0;
  prolongation[6](22,9) = -3.0/64.0;
  prolongation[6](22,11) = 1.0/64.0;
  prolongation[6](22,13) = 9.0/64.0;
  prolongation[6](22,15) = -3.0/64.0;
  prolongation[6](22,20) = -3.0/32.0;
  prolongation[6](22,21) = 9.0/32.0;
  prolongation[6](22,23) = 9.0/32.0;
  prolongation[6](22,25) = -3.0/32.0;
  prolongation[6](22,26) = 9.0/16.0;
  prolongation[6](23,1) = 1.0/64.0;
  prolongation[6](23,2) = -3.0/64.0;
  prolongation[6](23,5) = -3.0/64.0;
  prolongation[6](23,6) = 9.0/64.0;
  prolongation[6](23,9) = -3.0/32.0;
  prolongation[6](23,13) = 9.0/32.0;
  prolongation[6](23,17) = -3.0/32.0;
  prolongation[6](23,18) = 9.0/32.0;
  prolongation[6](23,23) = 9.0/16.0;
  prolongation[6](24,2) = -3.0/64.0;
  prolongation[6](24,3) = 1.0/64.0;
  prolongation[6](24,6) = 9.0/64.0;
  prolongation[6](24,7) = -3.0/64.0;
  prolongation[6](24,10) = -3.0/32.0;
  prolongation[6](24,14) = 9.0/32.0;
  prolongation[6](24,18) = 9.0/32.0;
  prolongation[6](24,19) = -3.0/32.0;
  prolongation[6](24,24) = 9.0/16.0;
  prolongation[6](25,8) = 1.0/64.0;
  prolongation[6](25,10) = -3.0/64.0;
  prolongation[6](25,12) = -3.0/64.0;
  prolongation[6](25,14) = 9.0/64.0;
  prolongation[6](25,20) = -3.0/32.0;
  prolongation[6](25,21) = 9.0/32.0;
  prolongation[6](25,22) = -3.0/32.0;
  prolongation[6](25,24) = 9.0/32.0;
  prolongation[6](25,26) = 9.0/16.0;
  prolongation[6](26,0) = -1.0/512.0;
  prolongation[6](26,1) = 3.0/512.0;
  prolongation[6](26,2) = -9.0/512.0;
  prolongation[6](26,3) = 3.0/512.0;
  prolongation[6](26,4) = 3.0/512.0;
  prolongation[6](26,5) = -9.0/512.0;
  prolongation[6](26,6) = 27.0/512.0;
  prolongation[6](26,7) = -9.0/512.0;
  prolongation[6](26,8) = 3.0/256.0;
  prolongation[6](26,9) = -9.0/256.0;
  prolongation[6](26,10) = -9.0/256.0;
  prolongation[6](26,11) = 3.0/256.0;
  prolongation[6](26,12) = -9.0/256.0;
  prolongation[6](26,13) = 27.0/256.0;
  prolongation[6](26,14) = 27.0/256.0;
  prolongation[6](26,15) = -9.0/256.0;
  prolongation[6](26,16) = 3.0/256.0;
  prolongation[6](26,17) = -9.0/256.0;
  prolongation[6](26,18) = 27.0/256.0;
  prolongation[6](26,19) = -9.0/256.0;
  prolongation[6](26,20) = -9.0/128.0;
  prolongation[6](26,21) = 27.0/128.0;
  prolongation[6](26,22) = -9.0/128.0;
  prolongation[6](26,23) = 27.0/128.0;
  prolongation[6](26,24) = 27.0/128.0;
  prolongation[6](26,25) = -9.0/128.0;
  prolongation[6](26,26) = 27.0/64.0;
  prolongation[7](0,25) = 1.0;
  prolongation[7](1,26) = 1.0;
  prolongation[7](2,24) = 1.0;
  prolongation[7](3,19) = 1.0;
  prolongation[7](4,15) = 1.0;
  prolongation[7](5,21) = 1.0;
  prolongation[7](6,14) = 1.0;
  prolongation[7](7,7) = 1.0;
  prolongation[7](8,23) = -1.0/8.0;
  prolongation[7](8,25) = 3.0/8.0;
  prolongation[7](8,26) = 3.0/4.0;
  prolongation[7](9,22) = -1.0/8.0;
  prolongation[7](9,24) = 3.0/8.0;
  prolongation[7](9,26) = 3.0/4.0;
  prolongation[7](10,18) = -1.0/8.0;
  prolongation[7](10,19) = 3.0/8.0;
  prolongation[7](10,24) = 3.0/4.0;
  prolongation[7](11,16) = -1.0/8.0;
  prolongation[7](11,19) = 3.0/8.0;
  prolongation[7](11,25) = 3.0/4.0;
  prolongation[7](12,13) = -1.0/8.0;
  prolongation[7](12,15) = 3.0/8.0;
  prolongation[7](12,21) = 3.0/4.0;
  prolongation[7](13,12) = -1.0/8.0;
  prolongation[7](13,14) = 3.0/8.0;
  prolongation[7](13,21) = 3.0/4.0;
  prolongation[7](14,6) = -1.0/8.0;
  prolongation[7](14,7) = 3.0/8.0;
  prolongation[7](14,14) = 3.0/4.0;
  prolongation[7](15,4) = -1.0/8.0;
  prolongation[7](15,7) = 3.0/8.0;
  prolongation[7](15,15) = 3.0/4.0;
  prolongation[7](16,11) = -1.0/8.0;
  prolongation[7](16,15) = 3.0/8.0;
  prolongation[7](16,25) = 3.0/4.0;
  prolongation[7](17,20) = -1.0/8.0;
  prolongation[7](17,21) = 3.0/8.0;
  prolongation[7](17,26) = 3.0/4.0;
  prolongation[7](18,10) = -1.0/8.0;
  prolongation[7](18,14) = 3.0/8.0;
  prolongation[7](18,24) = 3.0/4.0;
  prolongation[7](19,3) = -1.0/8.0;
  prolongation[7](19,7) = 3.0/8.0;
  prolongation[7](19,19) = 3.0/4.0;
  prolongation[7](20,16) = -3.0/64.0;
  prolongation[7](20,17) = 1.0/64.0;
  prolongation[7](20,18) = -3.0/64.0;
  prolongation[7](20,19) = 9.0/64.0;
  prolongation[7](20,22) = -3.0/32.0;
  prolongation[7](20,23) = -3.0/32.0;
  prolongation[7](20,24) = 9.0/32.0;
  prolongation[7](20,25) = 9.0/32.0;
  prolongation[7](20,26) = 9.0/16.0;
  prolongation[7](21,4) = -3.0/64.0;
  prolongation[7](21,5) = 1.0/64.0;
  prolongation[7](21,6) = -3.0/64.0;
  prolongation[7](21,7) = 9.0/64.0;
  prolongation[7](21,12) = -3.0/32.0;
  prolongation[7](21,13) = -3.0/32.0;
  prolongation[7](21,14) = 9.0/32.0;
  prolongation[7](21,15) = 9.0/32.0;
  prolongation[7](21,21) = 9.0/16.0;
  prolongation[7](22,9) = 1.0/64.0;
  prolongation[7](22,11) = -3.0/64.0;
  prolongation[7](22,13) = -3.0/64.0;
  prolongation[7](22,15) = 9.0/64.0;
  prolongation[7](22,20) = -3.0/32.0;
  prolongation[7](22,21) = 9.0/32.0;
  prolongation[7](22,23) = -3.0/32.0;
  prolongation[7](22,25) = 9.0/32.0;
  prolongation[7](22,26) = 9.0/16.0;
  prolongation[7](23,8) = 1.0/64.0;
  prolongation[7](23,10) = -3.0/64.0;
  prolongation[7](23,12) = -3.0/64.0;
  prolongation[7](23,14) = 9.0/64.0;
  prolongation[7](23,20) = -3.0/32.0;
  prolongation[7](23,21) = 9.0/32.0;
  prolongation[7](23,22) = -3.0/32.0;
  prolongation[7](23,24) = 9.0/32.0;
  prolongation[7](23,26) = 9.0/16.0;
  prolongation[7](24,2) = 1.0/64.0;
  prolongation[7](24,3) = -3.0/64.0;
  prolongation[7](24,6) = -3.0/64.0;
  prolongation[7](24,7) = 9.0/64.0;
  prolongation[7](24,10) = -3.0/32.0;
  prolongation[7](24,14) = 9.0/32.0;
  prolongation[7](24,18) = -3.0/32.0;
  prolongation[7](24,19) = 9.0/32.0;
  prolongation[7](24,24) = 9.0/16.0;
  prolongation[7](25,0) = 1.0/64.0;
  prolongation[7](25,3) = -3.0/64.0;
  prolongation[7](25,4) = -3.0/64.0;
  prolongation[7](25,7) = 9.0/64.0;
  prolongation[7](25,11) = -3.0/32.0;
  prolongation[7](25,15) = 9.0/32.0;
  prolongation[7](25,16) = -3.0/32.0;
  prolongation[7](25,19) = 9.0/32.0;
  prolongation[7](25,25) = 9.0/16.0;
  prolongation[7](26,0) = 3.0/512.0;
  prolongation[7](26,1) = -1.0/512.0;
  prolongation[7](26,2) = 3.0/512.0;
  prolongation[7](26,3) = -9.0/512.0;
  prolongation[7](26,4) = -9.0/512.0;
  prolongation[7](26,5) = 3.0/512.0;
  prolongation[7](26,6) = -9.0/512.0;
  prolongation[7](26,7) = 27.0/512.0;
  prolongation[7](26,8) = 3.0/256.0;
  prolongation[7](26,9) = 3.0/256.0;
  prolongation[7](26,10) = -9.0/256.0;
  prolongation[7](26,11) = -9.0/256.0;
  prolongation[7](26,12) = -9.0/256.0;
  prolongation[7](26,13) = -9.0/256.0;
  prolongation[7](26,14) = 27.0/256.0;
  prolongation[7](26,15) = 27.0/256.0;
  prolongation[7](26,16) = -9.0/256.0;
  prolongation[7](26,17) = 3.0/256.0;
  prolongation[7](26,18) = -9.0/256.0;
  prolongation[7](26,19) = 27.0/256.0;
  prolongation[7](26,20) = -9.0/128.0;
  prolongation[7](26,21) = 27.0/128.0;
  prolongation[7](26,22) = -9.0/128.0;
  prolongation[7](26,23) = -9.0/128.0;
  prolongation[7](26,24) = 27.0/128.0;
  prolongation[7](26,25) = 27.0/128.0;
  prolongation[7](26,26) = 27.0/64.0;

  restriction[0](0,0) = 1.0;
  restriction[0](8,1) = 1.0;
  restriction[0](11,3) = 1.0;
  restriction[0](16,4) = 1.0;
  restriction[0](20,2) = 1.0;
  restriction[0](22,5) = 1.0;
  restriction[0](25,7) = 1.0;
  restriction[0](26,6) = 1.0;
  restriction[1](1,1) = 1.0;
  restriction[1](8,0) = 1.0;
  restriction[1](9,2) = 1.0;
  restriction[1](17,5) = 1.0;
  restriction[1](20,3) = 1.0;
  restriction[1](22,4) = 1.0;
  restriction[1](23,6) = 1.0;
  restriction[1](26,7) = 1.0;
  restriction[2](2,2) = 1.0;
  restriction[2](9,1) = 1.0;
  restriction[2](10,3) = 1.0;
  restriction[2](18,6) = 1.0;
  restriction[2](20,0) = 1.0;
  restriction[2](23,5) = 1.0;
  restriction[2](24,7) = 1.0;
  restriction[2](26,4) = 1.0;
  restriction[3](3,3) = 1.0;
  restriction[3](10,2) = 1.0;
  restriction[3](11,0) = 1.0;
  restriction[3](19,7) = 1.0;
  restriction[3](20,1) = 1.0;
  restriction[3](24,6) = 1.0;
  restriction[3](25,4) = 1.0;
  restriction[3](26,5) = 1.0;
  restriction[4](4,4) = 1.0;
  restriction[4](12,5) = 1.0;
  restriction[4](15,7) = 1.0;
  restriction[4](16,0) = 1.0;
  restriction[4](21,6) = 1.0;
  restriction[4](22,1) = 1.0;
  restriction[4](25,3) = 1.0;
  restriction[4](26,2) = 1.0;
  restriction[5](5,5) = 1.0;
  restriction[5](12,4) = 1.0;
  restriction[5](13,6) = 1.0;
  restriction[5](17,1) = 1.0;
  restriction[5](21,7) = 1.0;
  restriction[5](22,0) = 1.0;
  restriction[5](23,2) = 1.0;
  restriction[5](26,3) = 1.0;
  restriction[6](6,6) = 1.0;
  restriction[6](13,5) = 1.0;
  restriction[6](14,7) = 1.0;
  restriction[6](18,2) = 1.0;
  restriction[6](21,4) = 1.0;
  restriction[6](23,1) = 1.0;
  restriction[6](24,3) = 1.0;
  restriction[6](26,0) = 1.0;
  restriction[7](7,7) = 1.0;
  restriction[7](14,6) = 1.0;
  restriction[7](15,4) = 1.0;
  restriction[7](19,3) = 1.0;
  restriction[7](21,5) = 1.0;
  restriction[7](24,2) = 1.0;
  restriction[7](25,0) = 1.0;
  restriction[7](26,1) = 1.0;
};



template <>
double
FEQuadraticSub<3>::shape_value (const unsigned int i,
				const Point<3>    &p) const
{
  Assert (i<total_dofs, ExcInvalidIndex(i));

  const double xi  = p(0),
	       eta = p(1),
	       zeta=p(2);
  switch (i)
    {
      case 0: return 1.0-3.0*xi+2.0*xi*xi+(-3.0+9.0*xi-6.0*xi*xi)*eta+(2.0
-6.0*xi+4.0*xi*xi)*eta*eta+(-3.0+9.0*xi-6.0*xi*xi+(9.0-27.0*xi+18.0*xi*xi)*eta+
(-6.0+18.0*xi-12.0*xi*xi)*eta*eta)*zeta+(2.0-6.0*xi+4.0*xi*xi+(-6.0+18.0*xi
-12.0*xi*xi)*eta+(4.0-12.0*xi+8.0*xi*xi)*eta*eta)*zeta*zeta;
      case 1: return -xi+2.0*xi*xi+(3.0*xi-6.0*xi*xi)*eta+(-2.0*xi+4.0*xi*xi)
*eta*eta+(3.0*xi-6.0*xi*xi+(-9.0*xi+18.0*xi*xi)*eta+(6.0*xi-12.0*xi*xi)*eta*eta
)*zeta+(-2.0*xi+4.0*xi*xi+(6.0*xi-12.0*xi*xi)*eta+(-4.0*xi+8.0*xi*xi)*eta*eta)*
zeta*zeta;
      case 2: return (xi-2.0*xi*xi+(-3.0*xi+6.0*xi*xi)*eta+(2.0*xi-4.0*xi*xi)
*eta*eta)*zeta+(-2.0*xi+4.0*xi*xi+(6.0*xi-12.0*xi*xi)*eta+(-4.0*xi+8.0*xi*xi)*
eta*eta)*zeta*zeta;
      case 3: return (-1.0+3.0*xi-2.0*xi*xi+(3.0-9.0*xi+6.0*xi*xi)*eta+(-2.0+
6.0*xi-4.0*xi*xi)*eta*eta)*zeta+(2.0-6.0*xi+4.0*xi*xi+(-6.0+18.0*xi-12.0*xi*xi)
*eta+(4.0-12.0*xi+8.0*xi*xi)*eta*eta)*zeta*zeta;
      case 4: return (-1.0+3.0*xi-2.0*xi*xi)*eta+(2.0-6.0*xi+4.0*xi*xi)*eta*
eta+((3.0-9.0*xi+6.0*xi*xi)*eta+(-6.0+18.0*xi-12.0*xi*xi)*eta*eta)*zeta+((-2.0+
6.0*xi-4.0*xi*xi)*eta+(4.0-12.0*xi+8.0*xi*xi)*eta*eta)*zeta*zeta;
      case 5: return (xi-2.0*xi*xi)*eta+(-2.0*xi+4.0*xi*xi)*eta*eta+((-3.0*xi
+6.0*xi*xi)*eta+(6.0*xi-12.0*xi*xi)*eta*eta)*zeta+((2.0*xi-4.0*xi*xi)*eta+(-4.0
*xi+8.0*xi*xi)*eta*eta)*zeta*zeta;
      case 6: return ((-xi+2.0*xi*xi)*eta+(2.0*xi-4.0*xi*xi)*eta*eta)*zeta+((
2.0*xi-4.0*xi*xi)*eta+(-4.0*xi+8.0*xi*xi)*eta*eta)*zeta*zeta;
      case 7: return ((1.0-3.0*xi+2.0*xi*xi)*eta+(-2.0+6.0*xi-4.0*xi*xi)*eta*
eta)*zeta+((-2.0+6.0*xi-4.0*xi*xi)*eta+(4.0-12.0*xi+8.0*xi*xi)*eta*eta)*zeta*
zeta;
      case 8: return 4.0*xi-4.0*xi*xi+(-12.0*xi+12.0*xi*xi)*eta+(8.0*xi-8.0*
xi*xi)*eta*eta+(-12.0*xi+12.0*xi*xi+(36.0*xi-36.0*xi*xi)*eta+(-24.0*xi+24.0*xi*
xi)*eta*eta)*zeta+(8.0*xi-8.0*xi*xi+(-24.0*xi+24.0*xi*xi)*eta+(16.0*xi-16.0*xi*
xi)*eta*eta)*zeta*zeta;
      case 9: return (-4.0*xi+8.0*xi*xi+(12.0*xi-24.0*xi*xi)*eta+(-8.0*xi+
16.0*xi*xi)*eta*eta)*zeta+(4.0*xi-8.0*xi*xi+(-12.0*xi+24.0*xi*xi)*eta+(8.0*xi
-16.0*xi*xi)*eta*eta)*zeta*zeta;
      case 10: return (-4.0*xi+4.0*xi*xi+(12.0*xi-12.0*xi*xi)*eta+(-8.0*xi+
8.0*xi*xi)*eta*eta)*zeta+(8.0*xi-8.0*xi*xi+(-24.0*xi+24.0*xi*xi)*eta+(16.0*xi
-16.0*xi*xi)*eta*eta)*zeta*zeta;
      case 11: return (4.0-12.0*xi+8.0*xi*xi+(-12.0+36.0*xi-24.0*xi*xi)*eta+(
8.0-24.0*xi+16.0*xi*xi)*eta*eta)*zeta+(-4.0+12.0*xi-8.0*xi*xi+(12.0-36.0*xi+
24.0*xi*xi)*eta+(-8.0+24.0*xi-16.0*xi*xi)*eta*eta)*zeta*zeta;
      case 12: return (-4.0*xi+4.0*xi*xi)*eta+(8.0*xi-8.0*xi*xi)*eta*eta+((
12.0*xi-12.0*xi*xi)*eta+(-24.0*xi+24.0*xi*xi)*eta*eta)*zeta+((-8.0*xi+8.0*xi*xi
)*eta+(16.0*xi-16.0*xi*xi)*eta*eta)*zeta*zeta;
      case 13: return ((4.0*xi-8.0*xi*xi)*eta+(-8.0*xi+16.0*xi*xi)*eta*eta)*
zeta+((-4.0*xi+8.0*xi*xi)*eta+(8.0*xi-16.0*xi*xi)*eta*eta)*zeta*zeta;
      case 14: return ((4.0*xi-4.0*xi*xi)*eta+(-8.0*xi+8.0*xi*xi)*eta*eta)*
zeta+((-8.0*xi+8.0*xi*xi)*eta+(16.0*xi-16.0*xi*xi)*eta*eta)*zeta*zeta;
      case 15: return ((-4.0+12.0*xi-8.0*xi*xi)*eta+(8.0-24.0*xi+16.0*xi*xi)*
eta*eta)*zeta+((4.0-12.0*xi+8.0*xi*xi)*eta+(-8.0+24.0*xi-16.0*xi*xi)*eta*eta)*
zeta*zeta;
      case 16: return (4.0-12.0*xi+8.0*xi*xi)*eta+(-4.0+12.0*xi-8.0*xi*xi)*
eta*eta+((-12.0+36.0*xi-24.0*xi*xi)*eta+(12.0-36.0*xi+24.0*xi*xi)*eta*eta)*zeta
+((8.0-24.0*xi+16.0*xi*xi)*eta+(-8.0+24.0*xi-16.0*xi*xi)*eta*eta)*zeta*zeta;
      case 17: return (-4.0*xi+8.0*xi*xi)*eta+(4.0*xi-8.0*xi*xi)*eta*eta+((
12.0*xi-24.0*xi*xi)*eta+(-12.0*xi+24.0*xi*xi)*eta*eta)*zeta+((-8.0*xi+16.0*xi*
xi)*eta+(8.0*xi-16.0*xi*xi)*eta*eta)*zeta*zeta;
      case 18: return ((4.0*xi-8.0*xi*xi)*eta+(-4.0*xi+8.0*xi*xi)*eta*eta)*
zeta+((-8.0*xi+16.0*xi*xi)*eta+(8.0*xi-16.0*xi*xi)*eta*eta)*zeta*zeta;
      case 19: return ((-4.0+12.0*xi-8.0*xi*xi)*eta+(4.0-12.0*xi+8.0*xi*xi)*
eta*eta)*zeta+((8.0-24.0*xi+16.0*xi*xi)*eta+(-8.0+24.0*xi-16.0*xi*xi)*eta*eta)*
zeta*zeta;
      case 20: return (16.0*xi-16.0*xi*xi+(-48.0*xi+48.0*xi*xi)*eta+(32.0*xi
-32.0*xi*xi)*eta*eta)*zeta+(-16.0*xi+16.0*xi*xi+(48.0*xi-48.0*xi*xi)*eta+(-32.0
*xi+32.0*xi*xi)*eta*eta)*zeta*zeta;
      case 21: return ((-16.0*xi+16.0*xi*xi)*eta+(32.0*xi-32.0*xi*xi)*eta*eta
)*zeta+((16.0*xi-16.0*xi*xi)*eta+(-32.0*xi+32.0*xi*xi)*eta*eta)*zeta*zeta;
      case 22: return (16.0*xi-16.0*xi*xi)*eta+(-16.0*xi+16.0*xi*xi)*eta*eta+
((-48.0*xi+48.0*xi*xi)*eta+(48.0*xi-48.0*xi*xi)*eta*eta)*zeta+((32.0*xi-32.0*xi
*xi)*eta+(-32.0*xi+32.0*xi*xi)*eta*eta)*zeta*zeta;
      case 23: return ((-16.0*xi+32.0*xi*xi)*eta+(16.0*xi-32.0*xi*xi)*eta*eta
)*zeta+((16.0*xi-32.0*xi*xi)*eta+(-16.0*xi+32.0*xi*xi)*eta*eta)*zeta*zeta;
      case 24: return ((-16.0*xi+16.0*xi*xi)*eta+(16.0*xi-16.0*xi*xi)*eta*eta
)*zeta+((32.0*xi-32.0*xi*xi)*eta+(-32.0*xi+32.0*xi*xi)*eta*eta)*zeta*zeta;
      case 25: return ((16.0-48.0*xi+32.0*xi*xi)*eta+(-16.0+48.0*xi-32.0*xi*
xi)*eta*eta)*zeta+((-16.0+48.0*xi-32.0*xi*xi)*eta+(16.0-48.0*xi+32.0*xi*xi)*eta
*eta)*zeta*zeta;
      case 26: return ((64.0*xi-64.0*xi*xi)*eta+(-64.0*xi+64.0*xi*xi)*eta*eta
)*zeta+((-64.0*xi+64.0*xi*xi)*eta+(64.0*xi-64.0*xi*xi)*eta*eta)*zeta*zeta;
    };
  return 0;
};



template <>
Tensor<1,3>
FEQuadraticSub<3>::shape_grad (const unsigned int i,
			       const Point<3>    &p) const
{
  Assert (i<total_dofs, ExcInvalidIndex(i));

  const double xi  = p(0),
	       eta = p(1),
	       zeta=p(2);
				   // originally, the return type of the
				   // function was Point<dim>, so we
				   // still construct it as that. it should
				   // make no difference in practice,
				   // however
  switch (i)
    {
      case 0: return Point<3>(-3.0+4.0*xi+(9.0-12.0*xi)*eta+(-6.0+8.0*xi)*eta*eta+(9.0-12.0*xi+(-27.0+36.0*xi)*eta+(18.0-24.0*xi)*eta*eta)*zeta+(-6.0+8.0*xi+(18.0-24.0*xi)*eta+(-12.0+16.0*xi)*eta*eta)*zeta*zeta,
      -3.0+9.0*xi-6.0*xi*xi+2.0*(2.0-6.0*xi+4.0*xi*xi)*eta+(9.0-27.0*xi+18.0*xi*xi+2.0*(-6.0+18.0*xi-12.0*xi*xi)*eta)*zeta+(-6.0+18.0*xi-12.0*xi*xi+2.0*(4.0-12.0*xi+8.0*xi*xi)*eta)*zeta*zeta,
      -3.0+9.0*xi-6.0*xi*xi+(9.0-27.0*xi+18.0*xi*xi)*eta+(-6.0+18.0*xi-12.0*xi*xi)*eta*eta+2.0*(2.0-6.0*xi+4.0*xi*xi+(-6.0+18.0*xi-12.0*xi*xi)*eta+(4.0-12.0*xi+8.0*xi*xi)*eta*eta)*zeta);
      case 1: return Point<3>(-1.0+4.0*xi+(3.0-12.0*xi)*eta+(-2.0+8.0*xi)*eta*eta+(3.0-12.0*xi+(-9.0+36.0*xi)*eta+(6.0-24.0*xi)*eta*eta)*zeta+(-2.0+8.0*xi+(6.0-24.0*xi)*eta+(-4.0+16.0*xi)*eta*eta)*zeta*zeta,
      3.0*xi-6.0*xi*xi+2.0*(-2.0*xi+4.0*xi*xi)*eta+(-9.0*xi+18.0*xi*xi+2.0*(6.0*xi-12.0*xi*xi)*eta)*zeta+(6.0*xi-12.0*xi*xi+2.0*(-4.0*xi+8.0*xi*xi)*eta)*zeta*zeta,
      3.0*xi-6.0*xi*xi+(-9.0*xi+18.0*xi*xi)*eta+(6.0*xi-12.0*xi*xi)*eta*eta+2.0*(-2.0*xi+4.0*xi*xi+(6.0*xi-12.0*xi*xi)*eta+(-4.0*xi+8.0*xi*xi)*eta*eta)*zeta);
      case 2: return Point<3>((1.0-4.0*xi+(-3.0+12.0*xi)*eta+(2.0-8.0*xi)*eta*eta)*zeta+(-2.0+8.0*xi+(6.0-24.0*xi)*eta+(-4.0+16.0*xi)*eta*eta)*zeta*zeta,
      (-3.0*xi+6.0*xi*xi+2.0*(2.0*xi-4.0*xi*xi)*eta)*zeta+(6.0*xi-12.0*xi*xi+2.0*(-4.0*xi+8.0*xi*xi)*eta)*zeta*zeta,
      xi-2.0*xi*xi+(-3.0*xi+6.0*xi*xi)*eta+(2.0*xi-4.0*xi*xi)*eta*eta+2.0*(-2.0*xi+4.0*xi*xi+(6.0*xi-12.0*xi*xi)*eta+(-4.0*xi+8.0*xi*xi)*eta*eta)*zeta);
      case 3: return Point<3>((3.0-4.0*xi+(-9.0+12.0*xi)*eta+(6.0-8.0*xi)*eta*eta)*zeta+(-6.0+8.0*xi+(18.0-24.0*xi)*eta+(-12.0+16.0*xi)*eta*eta)*zeta*zeta,
      (3.0-9.0*xi+6.0*xi*xi+2.0*(-2.0+6.0*xi-4.0*xi*xi)*eta)*zeta+(-6.0+18.0*xi-12.0*xi*xi+2.0*(4.0-12.0*xi+8.0*xi*xi)*eta)*zeta*zeta,
      -1.0+3.0*xi-2.0*xi*xi+(3.0-9.0*xi+6.0*xi*xi)*eta+(-2.0+6.0*xi-4.0*xi*xi)*eta*eta+2.0*(2.0-6.0*xi+4.0*xi*xi+(-6.0+18.0*xi-12.0*xi*xi)*eta+(4.0-12.0*xi+8.0*xi*xi)*eta*eta)*zeta);
      case 4: return Point<3>((3.0-4.0*xi)*eta+(-6.0+8.0*xi)*eta*eta+((-9.0+12.0*xi)*eta+(18.0-24.0*xi)*eta*eta)*zeta+((6.0-8.0*xi)*eta+(-12.0+16.0*xi)*eta*eta)*zeta*zeta,
      -1.0+3.0*xi-2.0*xi*xi+2.0*(2.0-6.0*xi+4.0*xi*xi)*eta+(3.0-9.0*xi+6.0*xi*xi+2.0*(-6.0+18.0*xi-12.0*xi*xi)*eta)*zeta+(-2.0+6.0*xi-4.0*xi*xi+2.0*(4.0-12.0*xi+8.0*xi*xi)*eta)*zeta*zeta,
      (3.0-9.0*xi+6.0*xi*xi)*eta+(-6.0+18.0*xi-12.0*xi*xi)*eta*eta+2.0*((-2.0+6.0*xi-4.0*xi*xi)*eta+(4.0-12.0*xi+8.0*xi*xi)*eta*eta)*zeta);
      case 5: return Point<3>((1.0-4.0*xi)*eta+(-2.0+8.0*xi)*eta*eta+((-3.0+12.0*xi)*eta+(6.0-24.0*xi)*eta*eta)*zeta+((2.0-8.0*xi)*eta+(-4.0+16.0*xi)*eta*eta)*zeta*zeta,
      xi-2.0*xi*xi+2.0*(-2.0*xi+4.0*xi*xi)*eta+(-3.0*xi+6.0*xi*xi+2.0*(6.0*xi-12.0*xi*xi)*eta)*zeta+(2.0*xi-4.0*xi*xi+2.0*(-4.0*xi+8.0*xi*xi)*eta)*zeta*zeta,
      (-3.0*xi+6.0*xi*xi)*eta+(6.0*xi-12.0*xi*xi)*eta*eta+2.0*((2.0*xi-4.0*xi*xi)*eta+(-4.0*xi+8.0*xi*xi)*eta*eta)*zeta);
      case 6: return Point<3>(((-1.0+4.0*xi)*eta+(2.0-8.0*xi)*eta*eta)*zeta+((2.0-8.0*xi)*eta+(-4.0+16.0*xi)*eta*eta)*zeta*zeta,
      (-xi+2.0*xi*xi+2.0*(2.0*xi-4.0*xi*xi)*eta)*zeta+(2.0*xi-4.0*xi*xi+2.0*(-4.0*xi+8.0*xi*xi)*eta)*zeta*zeta,
      (-xi+2.0*xi*xi)*eta+(2.0*xi-4.0*xi*xi)*eta*eta+2.0*((2.0*xi-4.0*xi*xi)*eta+(-4.0*xi+8.0*xi*xi)*eta*eta)*zeta);
      case 7: return Point<3>(((-3.0+4.0*xi)*eta+(6.0-8.0*xi)*eta*eta)*zeta+((6.0-8.0*xi)*eta+(-12.0+16.0*xi)*eta*eta)*zeta*zeta,
      (1.0-3.0*xi+2.0*xi*xi+2.0*(-2.0+6.0*xi-4.0*xi*xi)*eta)*zeta+(-2.0+6.0*xi-4.0*xi*xi+2.0*(4.0-12.0*xi+8.0*xi*xi)*eta)*zeta*zeta,
      (1.0-3.0*xi+2.0*xi*xi)*eta+(-2.0+6.0*xi-4.0*xi*xi)*eta*eta+2.0*((-2.0+6.0*xi-4.0*xi*xi)*eta+(4.0-12.0*xi+8.0*xi*xi)*eta*eta)*zeta);
      case 8: return Point<3>(4.0-8.0*xi+(-12.0+24.0*xi)*eta+(8.0-16.0*xi)*eta*eta+(-12.0+24.0*xi+(36.0-72.0*xi)*eta+(-24.0+48.0*xi)*eta*eta)*zeta+(8.0-16.0*xi+(-24.0+48.0*xi)*eta+(16.0-32.0*xi)*eta*eta)*zeta*zeta,
      -12.0*xi+12.0*xi*xi+2.0*(8.0*xi-8.0*xi*xi)*eta+(36.0*xi-36.0*xi*xi+2.0*(-24.0*xi+24.0*xi*xi)*eta)*zeta+(-24.0*xi+24.0*xi*xi+2.0*(16.0*xi-16.0*xi*xi)*eta)*zeta*zeta,
      -12.0*xi+12.0*xi*xi+(36.0*xi-36.0*xi*xi)*eta+(-24.0*xi+24.0*xi*xi)*eta*eta+2.0*(8.0*xi-8.0*xi*xi+(-24.0*xi+24.0*xi*xi)*eta+(16.0*xi-16.0*xi*xi)*eta*eta)*zeta);
      case 9: return Point<3>((-4.0+16.0*xi+(12.0-48.0*xi)*eta+(-8.0+32.0*xi)*eta*eta)*zeta+(4.0-16.0*xi+(-12.0+48.0*xi)*eta+(8.0-32.0*xi)*eta*eta)*zeta*zeta,
      (12.0*xi-24.0*xi*xi+2.0*(-8.0*xi+16.0*xi*xi)*eta)*zeta+(-12.0*xi+24.0*xi*xi+2.0*(8.0*xi-16.0*xi*xi)*eta)*zeta*zeta,
      -4.0*xi+8.0*xi*xi+(12.0*xi-24.0*xi*xi)*eta+(-8.0*xi+16.0*xi*xi)*eta*eta+2.0*(4.0*xi-8.0*xi*xi+(-12.0*xi+24.0*xi*xi)*eta+(8.0*xi-16.0*xi*xi)*eta*eta)*zeta);
      case 10: return Point<3>((-4.0+8.0*xi+(12.0-24.0*xi)*eta+(-8.0+16.0*xi)*eta*eta)*zeta+(8.0-16.0*xi+(-24.0+48.0*xi)*eta+(16.0-32.0*xi)*eta*eta)*zeta*zeta,
      (12.0*xi-12.0*xi*xi+2.0*(-8.0*xi+8.0*xi*xi)*eta)*zeta+(-24.0*xi+24.0*xi*xi+2.0*(16.0*xi-16.0*xi*xi)*eta)*zeta*zeta,
      -4.0*xi+4.0*xi*xi+(12.0*xi-12.0*xi*xi)*eta+(-8.0*xi+8.0*xi*xi)*eta*eta+2.0*(8.0*xi-8.0*xi*xi+(-24.0*xi+24.0*xi*xi)*eta+(16.0*xi-16.0*xi*xi)*eta*eta)*zeta);
      case 11: return Point<3>((-12.0+16.0*xi+(36.0-48.0*xi)*eta+(-24.0+32.0*xi)*eta*eta)*zeta+(12.0-16.0*xi+(-36.0+48.0*xi)*eta+(24.0-32.0*xi)*eta*eta)*zeta*zeta,
      (-12.0+36.0*xi-24.0*xi*xi+2.0*(8.0-24.0*xi+16.0*xi*xi)*eta)*zeta+(12.0-36.0*xi+24.0*xi*xi+2.0*(-8.0+24.0*xi-16.0*xi*xi)*eta)*zeta*zeta,
      4.0-12.0*xi+8.0*xi*xi+(-12.0+36.0*xi-24.0*xi*xi)*eta+(8.0-24.0*xi+16.0*xi*xi)*eta*eta+2.0*(-4.0+12.0*xi-8.0*xi*xi+(12.0-36.0*xi+24.0*xi*xi)*eta+(-8.0+24.0*xi-16.0*xi*xi)*eta*eta)*zeta);
      case 12: return Point<3>((-4.0+8.0*xi)*eta+(8.0-16.0*xi)*eta*eta+((12.0-24.0*xi)*eta+(-24.0+48.0*xi)*eta*eta)*zeta+((-8.0+16.0*xi)*eta+(16.0-32.0*xi)*eta*eta)*zeta*zeta,
      -4.0*xi+4.0*xi*xi+2.0*(8.0*xi-8.0*xi*xi)*eta+(12.0*xi-12.0*xi*xi+2.0*(-24.0*xi+24.0*xi*xi)*eta)*zeta+(-8.0*xi+8.0*xi*xi+2.0*(16.0*xi-16.0*xi*xi)*eta)*zeta*zeta,
      (12.0*xi-12.0*xi*xi)*eta+(-24.0*xi+24.0*xi*xi)*eta*eta+2.0*((-8.0*xi+8.0*xi*xi)*eta+(16.0*xi-16.0*xi*xi)*eta*eta)*zeta);
      case 13: return Point<3>(((4.0-16.0*xi)*eta+(-8.0+32.0*xi)*eta*eta)*zeta+((-4.0+16.0*xi)*eta+(8.0-32.0*xi)*eta*eta)*zeta*zeta,
      (4.0*xi-8.0*xi*xi+2.0*(-8.0*xi+16.0*xi*xi)*eta)*zeta+(-4.0*xi+8.0*xi*xi+2.0*(8.0*xi-16.0*xi*xi)*eta)*zeta*zeta,
      (4.0*xi-8.0*xi*xi)*eta+(-8.0*xi+16.0*xi*xi)*eta*eta+2.0*((-4.0*xi+8.0*xi*xi)*eta+(8.0*xi-16.0*xi*xi)*eta*eta)*zeta);
      case 14: return Point<3>(((4.0-8.0*xi)*eta+(-8.0+16.0*xi)*eta*eta)*zeta+((-8.0+16.0*xi)*eta+(16.0-32.0*xi)*eta*eta)*zeta*zeta,
      (4.0*xi-4.0*xi*xi+2.0*(-8.0*xi+8.0*xi*xi)*eta)*zeta+(-8.0*xi+8.0*xi*xi+2.0*(16.0*xi-16.0*xi*xi)*eta)*zeta*zeta,
      (4.0*xi-4.0*xi*xi)*eta+(-8.0*xi+8.0*xi*xi)*eta*eta+2.0*((-8.0*xi+8.0*xi*xi)*eta+(16.0*xi-16.0*xi*xi)*eta*eta)*zeta);
      case 15: return Point<3>(((12.0-16.0*xi)*eta+(-24.0+32.0*xi)*eta*eta)*zeta+((-12.0+16.0*xi)*eta+(24.0-32.0*xi)*eta*eta)*zeta*zeta,
      (-4.0+12.0*xi-8.0*xi*xi+2.0*(8.0-24.0*xi+16.0*xi*xi)*eta)*zeta+(4.0-12.0*xi+8.0*xi*xi+2.0*(-8.0+24.0*xi-16.0*xi*xi)*eta)*zeta*zeta,
      (-4.0+12.0*xi-8.0*xi*xi)*eta+(8.0-24.0*xi+16.0*xi*xi)*eta*eta+2.0*((4.0-12.0*xi+8.0*xi*xi)*eta+(-8.0+24.0*xi-16.0*xi*xi)*eta*eta)*zeta);
      case 16: return Point<3>((-12.0+16.0*xi)*eta+(12.0-16.0*xi)*eta*eta+((36.0-48.0*xi)*eta+(-36.0+48.0*xi)*eta*eta)*zeta+((-24.0+32.0*xi)*eta+(24.0-32.0*xi)*eta*eta)*zeta*zeta,
      4.0-12.0*xi+8.0*xi*xi+2.0*(-4.0+12.0*xi-8.0*xi*xi)*eta+(-12.0+36.0*xi-24.0*xi*xi+2.0*(12.0-36.0*xi+24.0*xi*xi)*eta)*zeta+(8.0-24.0*xi+16.0*xi*xi+2.0*(-8.0+24.0*xi-16.0*xi*xi)*eta)*zeta*zeta,
      (-12.0+36.0*xi-24.0*xi*xi)*eta+(12.0-36.0*xi+24.0*xi*xi)*eta*eta+2.0*((8.0-24.0*xi+16.0*xi*xi)*eta+(-8.0+24.0*xi-16.0*xi*xi)*eta*eta)*zeta);
      case 17: return Point<3>((-4.0+16.0*xi)*eta+(4.0-16.0*xi)*eta*eta+((12.0-48.0*xi)*eta+(-12.0+48.0*xi)*eta*eta)*zeta+((-8.0+32.0*xi)*eta+(8.0-32.0*xi)*eta*eta)*zeta*zeta,
      -4.0*xi+8.0*xi*xi+2.0*(4.0*xi-8.0*xi*xi)*eta+(12.0*xi-24.0*xi*xi+2.0*(-12.0*xi+24.0*xi*xi)*eta)*zeta+(-8.0*xi+16.0*xi*xi+2.0*(8.0*xi-16.0*xi*xi)*eta)*zeta*zeta,
      (12.0*xi-24.0*xi*xi)*eta+(-12.0*xi+24.0*xi*xi)*eta*eta+2.0*((-8.0*xi+16.0*xi*xi)*eta+(8.0*xi-16.0*xi*xi)*eta*eta)*zeta);
      case 18: return Point<3>(((4.0-16.0*xi)*eta+(-4.0+16.0*xi)*eta*eta)*zeta+((-8.0+32.0*xi)*eta+(8.0-32.0*xi)*eta*eta)*zeta*zeta,
      (4.0*xi-8.0*xi*xi+2.0*(-4.0*xi+8.0*xi*xi)*eta)*zeta+(-8.0*xi+16.0*xi*xi+2.0*(8.0*xi-16.0*xi*xi)*eta)*zeta*zeta,
      (4.0*xi-8.0*xi*xi)*eta+(-4.0*xi+8.0*xi*xi)*eta*eta+2.0*((-8.0*xi+16.0*xi*xi)*eta+(8.0*xi-16.0*xi*xi)*eta*eta)*zeta);
      case 19: return Point<3>(((12.0-16.0*xi)*eta+(-12.0+16.0*xi)*eta*eta)*zeta+((-24.0+32.0*xi)*eta+(24.0-32.0*xi)*eta*eta)*zeta*zeta,
      (-4.0+12.0*xi-8.0*xi*xi+2.0*(4.0-12.0*xi+8.0*xi*xi)*eta)*zeta+(8.0-24.0*xi+16.0*xi*xi+2.0*(-8.0+24.0*xi-16.0*xi*xi)*eta)*zeta*zeta,
      (-4.0+12.0*xi-8.0*xi*xi)*eta+(4.0-12.0*xi+8.0*xi*xi)*eta*eta+2.0*((8.0-24.0*xi+16.0*xi*xi)*eta+(-8.0+24.0*xi-16.0*xi*xi)*eta*eta)*zeta);
      case 20: return Point<3>((16.0-32.0*xi+(-48.0+96.0*xi)*eta+(32.0-64.0*xi)*eta*eta)*zeta+(-16.0+32.0*xi+(48.0-96.0*xi)*eta+(-32.0+64.0*xi)*eta*eta)*zeta*zeta,
      (-48.0*xi+48.0*xi*xi+2.0*(32.0*xi-32.0*xi*xi)*eta)*zeta+(48.0*xi-48.0*xi*xi+2.0*(-32.0*xi+32.0*xi*xi)*eta)*zeta*zeta,
      16.0*xi-16.0*xi*xi+(-48.0*xi+48.0*xi*xi)*eta+(32.0*xi-32.0*xi*xi)*eta*eta+2.0*(-16.0*xi+16.0*xi*xi+(48.0*xi-48.0*xi*xi)*eta+(-32.0*xi+32.0*xi*xi)*eta*eta)*zeta);
      case 21: return Point<3>(((-16.0+32.0*xi)*eta+(32.0-64.0*xi)*eta*eta)*zeta+((16.0-32.0*xi)*eta+(-32.0+64.0*xi)*eta*eta)*zeta*zeta,
      (-16.0*xi+16.0*xi*xi+2.0*(32.0*xi-32.0*xi*xi)*eta)*zeta+(16.0*xi-16.0*xi*xi+2.0*(-32.0*xi+32.0*xi*xi)*eta)*zeta*zeta,
      (-16.0*xi+16.0*xi*xi)*eta+(32.0*xi-32.0*xi*xi)*eta*eta+2.0*((16.0*xi-16.0*xi*xi)*eta+(-32.0*xi+32.0*xi*xi)*eta*eta)*zeta);
      case 22: return Point<3>((16.0-32.0*xi)*eta+(-16.0+32.0*xi)*eta*eta+((-48.0+96.0*xi)*eta+(48.0-96.0*xi)*eta*eta)*zeta+((32.0-64.0*xi)*eta+(-32.0+64.0*xi)*eta*eta)*zeta*zeta,
      16.0*xi-16.0*xi*xi+2.0*(-16.0*xi+16.0*xi*xi)*eta+(-48.0*xi+48.0*xi*xi+2.0*(48.0*xi-48.0*xi*xi)*eta)*zeta+(32.0*xi-32.0*xi*xi+2.0*(-32.0*xi+32.0*xi*xi)*eta)*zeta*zeta,
      (-48.0*xi+48.0*xi*xi)*eta+(48.0*xi-48.0*xi*xi)*eta*eta+2.0*((32.0*xi-32.0*xi*xi)*eta+(-32.0*xi+32.0*xi*xi)*eta*eta)*zeta);
      case 23: return Point<3>(((-16.0+64.0*xi)*eta+(16.0-64.0*xi)*eta*eta)*zeta+((16.0-64.0*xi)*eta+(-16.0+64.0*xi)*eta*eta)*zeta*zeta,
      (-16.0*xi+32.0*xi*xi+2.0*(16.0*xi-32.0*xi*xi)*eta)*zeta+(16.0*xi-32.0*xi*xi+2.0*(-16.0*xi+32.0*xi*xi)*eta)*zeta*zeta,
      (-16.0*xi+32.0*xi*xi)*eta+(16.0*xi-32.0*xi*xi)*eta*eta+2.0*((16.0*xi-32.0*xi*xi)*eta+(-16.0*xi+32.0*xi*xi)*eta*eta)*zeta);
      case 24: return Point<3>(((-16.0+32.0*xi)*eta+(16.0-32.0*xi)*eta*eta)*zeta+((32.0-64.0*xi)*eta+(-32.0+64.0*xi)*eta*eta)*zeta*zeta,
      (-16.0*xi+16.0*xi*xi+2.0*(16.0*xi-16.0*xi*xi)*eta)*zeta+(32.0*xi-32.0*xi*xi+2.0*(-32.0*xi+32.0*xi*xi)*eta)*zeta*zeta,
      (-16.0*xi+16.0*xi*xi)*eta+(16.0*xi-16.0*xi*xi)*eta*eta+2.0*((32.0*xi-32.0*xi*xi)*eta+(-32.0*xi+32.0*xi*xi)*eta*eta)*zeta);
      case 25: return Point<3>(((-48.0+64.0*xi)*eta+(48.0-64.0*xi)*eta*eta)*zeta+((48.0-64.0*xi)*eta+(-48.0+64.0*xi)*eta*eta)*zeta*zeta,
      (16.0-48.0*xi+32.0*xi*xi+2.0*(-16.0+48.0*xi-32.0*xi*xi)*eta)*zeta+(-16.0+48.0*xi-32.0*xi*xi+2.0*(16.0-48.0*xi+32.0*xi*xi)*eta)*zeta*zeta,
      (16.0-48.0*xi+32.0*xi*xi)*eta+(-16.0+48.0*xi-32.0*xi*xi)*eta*eta+2.0*((-16.0+48.0*xi-32.0*xi*xi)*eta+(16.0-48.0*xi+32.0*xi*xi)*eta*eta)*zeta);
      case 26: return Point<3>(((64.0-128.0*xi)*eta+(-64.0+128.0*xi)*eta*eta)*zeta+((-64.0+128.0*xi)*eta+(64.0-128.0*xi)*eta*eta)*zeta*zeta,
      (64.0*xi-64.0*xi*xi+2.0*(-64.0*xi+64.0*xi*xi)*eta)*zeta+(-64.0*xi+64.0*xi*xi+2.0*(64.0*xi-64.0*xi*xi)*eta)*zeta*zeta,
      (64.0*xi-64.0*xi*xi)*eta+(-64.0*xi+64.0*xi*xi)*eta*eta+2.0*((-64.0*xi+64.0*xi*xi)*eta+(64.0*xi-64.0*xi*xi)*eta*eta)*zeta);
    };
  return Point<3> ();
};



template <>
Tensor<2,3>
FEQuadraticSub<3>::shape_grad_grad (const unsigned int i,
				    const Point<3>    &p) const
{
  Assert (i<total_dofs, ExcInvalidIndex(i));

  const double xi  = p(0),
	       eta = p(1),
	       zeta=p(2);
  Tensor<2,3> return_value;
  
  switch (i)
    {
case 0:
return_value[0][0] = 4.0-12.0*eta+8.0*eta*eta+(-12.0+36.0*eta-24.0*eta*eta)*zeta+(8.0-24.0*eta+16.0*eta*eta)*zeta*zeta;
return_value[0][1] = 9.0-12.0*xi+2.0*(-6.0+8.0*xi)*eta+(-27.0+36.0*xi+2.0*(18.0-24.0*xi)*eta)*zeta+(18.0-24.0*xi+2.0*(-12.0+16.0*xi)*eta)*zeta*zeta;
return_value[0][2] = 9.0-12.0*xi+(-27.0+36.0*xi)*eta+(18.0-24.0*xi)*eta*eta+2.0*(-6.0+8.0*xi+(18.0-24.0*xi)*eta+(-12.0+16.0*xi)*eta*eta)*zeta;
return_value[1][0] = 9.0-12.0*xi+2.0*(-6.0+8.0*xi)*eta+(-27.0+36.0*xi+2.0*(18.0-24.0*xi)*eta)*zeta+(18.0-24.0*xi+2.0*(-12.0+16.0*xi)*eta)*zeta*zeta;
return_value[1][1] = 4.0-12.0*xi+8.0*xi*xi+(-12.0+36.0*xi-24.0*xi*xi)*zeta+(8.0-24.0*xi+16.0*xi*xi)*zeta*zeta;
return_value[1][2] = 9.0-27.0*xi+18.0*xi*xi+2.0*(-6.0+18.0*xi-12.0*xi*xi)*eta+2.0*(-6.0+18.0*xi-12.0*xi*xi+2.0*(4.0-12.0*xi+8.0*xi*xi)*eta)*zeta;
return_value[2][0] = 9.0-12.0*xi+(-27.0+36.0*xi)*eta+(18.0-24.0*xi)*eta*eta+2.0*(-6.0+8.0*xi+(18.0-24.0*xi)*eta+(-12.0+16.0*xi)*eta*eta)*zeta;
return_value[2][1] = 9.0-27.0*xi+18.0*xi*xi+2.0*(-6.0+18.0*xi-12.0*xi*xi)*eta+2.0*(-6.0+18.0*xi-12.0*xi*xi+2.0*(4.0-12.0*xi+8.0*xi*xi)*eta)*zeta;
return_value[2][2] = 4.0-12.0*xi+8.0*xi*xi+2.0*(-6.0+18.0*xi-12.0*xi*xi)*eta+2.0*(4.0-12.0*xi+8.0*xi*xi)*eta*eta;
break;
case 1:
return_value[0][0] = 4.0-12.0*eta+8.0*eta*eta+(-12.0+36.0*eta-24.0*eta*eta)*zeta+(8.0-24.0*eta+16.0*eta*eta)*zeta*zeta;
return_value[0][1] = 3.0-12.0*xi+2.0*(-2.0+8.0*xi)*eta+(-9.0+36.0*xi+2.0*(6.0-24.0*xi)*eta)*zeta+(6.0-24.0*xi+2.0*(-4.0+16.0*xi)*eta)*zeta*zeta;
return_value[0][2] = 3.0-12.0*xi+(-9.0+36.0*xi)*eta+(6.0-24.0*xi)*eta*eta+2.0*(-2.0+8.0*xi+(6.0-24.0*xi)*eta+(-4.0+16.0*xi)*eta*eta)*zeta;
return_value[1][0] = 3.0-12.0*xi+2.0*(-2.0+8.0*xi)*eta+(-9.0+36.0*xi+2.0*(6.0-24.0*xi)*eta)*zeta+(6.0-24.0*xi+2.0*(-4.0+16.0*xi)*eta)*zeta*zeta;
return_value[1][1] = -4.0*xi+8.0*xi*xi+(12.0*xi-24.0*xi*xi)*zeta+(-8.0*xi+16.0*xi*xi)*zeta*zeta;
return_value[1][2] = -9.0*xi+18.0*xi*xi+2.0*(6.0*xi-12.0*xi*xi)*eta+2.0*(6.0*xi-12.0*xi*xi+2.0*(-4.0*xi+8.0*xi*xi)*eta)*zeta;
return_value[2][0] = 3.0-12.0*xi+(-9.0+36.0*xi)*eta+(6.0-24.0*xi)*eta*eta+2.0*(-2.0+8.0*xi+(6.0-24.0*xi)*eta+(-4.0+16.0*xi)*eta*eta)*zeta;
return_value[2][1] = -9.0*xi+18.0*xi*xi+2.0*(6.0*xi-12.0*xi*xi)*eta+2.0*(6.0*xi-12.0*xi*xi+2.0*(-4.0*xi+8.0*xi*xi)*eta)*zeta;
return_value[2][2] = -4.0*xi+8.0*xi*xi+2.0*(6.0*xi-12.0*xi*xi)*eta+2.0*(-4.0*xi+8.0*xi*xi)*eta*eta;
break;
case 2:
return_value[0][0] = (-4.0+12.0*eta-8.0*eta*eta)*zeta+(8.0-24.0*eta+16.0*eta*eta)*zeta*zeta;
return_value[0][1] = (-3.0+12.0*xi+2.0*(2.0-8.0*xi)*eta)*zeta+(6.0-24.0*xi+2.0*(-4.0+16.0*xi)*eta)*zeta*zeta;
return_value[0][2] = 1.0-4.0*xi+(-3.0+12.0*xi)*eta+(2.0-8.0*xi)*eta*eta+2.0*(-2.0+8.0*xi+(6.0-24.0*xi)*eta+(-4.0+16.0*xi)*eta*eta)*zeta;
return_value[1][0] = (-3.0+12.0*xi+2.0*(2.0-8.0*xi)*eta)*zeta+(6.0-24.0*xi+2.0*(-4.0+16.0*xi)*eta)*zeta*zeta;
return_value[1][1] = (4.0*xi-8.0*xi*xi)*zeta+(-8.0*xi+16.0*xi*xi)*zeta*zeta;
return_value[1][2] = -3.0*xi+6.0*xi*xi+2.0*(2.0*xi-4.0*xi*xi)*eta+2.0*(6.0*xi-12.0*xi*xi+2.0*(-4.0*xi+8.0*xi*xi)*eta)*zeta;
return_value[2][0] = 1.0-4.0*xi+(-3.0+12.0*xi)*eta+(2.0-8.0*xi)*eta*eta+2.0*(-2.0+8.0*xi+(6.0-24.0*xi)*eta+(-4.0+16.0*xi)*eta*eta)*zeta;
return_value[2][1] = -3.0*xi+6.0*xi*xi+2.0*(2.0*xi-4.0*xi*xi)*eta+2.0*(6.0*xi-12.0*xi*xi+2.0*(-4.0*xi+8.0*xi*xi)*eta)*zeta;
return_value[2][2] = -4.0*xi+8.0*xi*xi+2.0*(6.0*xi-12.0*xi*xi)*eta+2.0*(-4.0*xi+8.0*xi*xi)*eta*eta;
break;
case 3:
return_value[0][0] = (-4.0+12.0*eta-8.0*eta*eta)*zeta+(8.0-24.0*eta+16.0*eta*eta)*zeta*zeta;
return_value[0][1] = (-9.0+12.0*xi+2.0*(6.0-8.0*xi)*eta)*zeta+(18.0-24.0*xi+2.0*(-12.0+16.0*xi)*eta)*zeta*zeta;
return_value[0][2] = 3.0-4.0*xi+(-9.0+12.0*xi)*eta+(6.0-8.0*xi)*eta*eta+2.0*(-6.0+8.0*xi+(18.0-24.0*xi)*eta+(-12.0+16.0*xi)*eta*eta)*zeta;
return_value[1][0] = (-9.0+12.0*xi+2.0*(6.0-8.0*xi)*eta)*zeta+(18.0-24.0*xi+2.0*(-12.0+16.0*xi)*eta)*zeta*zeta;
return_value[1][1] = (-4.0+12.0*xi-8.0*xi*xi)*zeta+(8.0-24.0*xi+16.0*xi*xi)*zeta*zeta;
return_value[1][2] = 3.0-9.0*xi+6.0*xi*xi+2.0*(-2.0+6.0*xi-4.0*xi*xi)*eta+2.0*(-6.0+18.0*xi-12.0*xi*xi+2.0*(4.0-12.0*xi+8.0*xi*xi)*eta)*zeta;
return_value[2][0] = 3.0-4.0*xi+(-9.0+12.0*xi)*eta+(6.0-8.0*xi)*eta*eta+2.0*(-6.0+8.0*xi+(18.0-24.0*xi)*eta+(-12.0+16.0*xi)*eta*eta)*zeta;
return_value[2][1] = 3.0-9.0*xi+6.0*xi*xi+2.0*(-2.0+6.0*xi-4.0*xi*xi)*eta+2.0*(-6.0+18.0*xi-12.0*xi*xi+2.0*(4.0-12.0*xi+8.0*xi*xi)*eta)*zeta;
return_value[2][2] = 4.0-12.0*xi+8.0*xi*xi+2.0*(-6.0+18.0*xi-12.0*xi*xi)*eta+2.0*(4.0-12.0*xi+8.0*xi*xi)*eta*eta;
break;
case 4:
return_value[0][0] = -4.0*eta+8.0*eta*eta+(12.0*eta-24.0*eta*eta)*zeta+(-8.0*eta+16.0*eta*eta)*zeta*zeta;
return_value[0][1] = 3.0-4.0*xi+2.0*(-6.0+8.0*xi)*eta+(-9.0+12.0*xi+2.0*(18.0-24.0*xi)*eta)*zeta+(6.0-8.0*xi+2.0*(-12.0+16.0*xi)*eta)*zeta*zeta;
return_value[0][2] = (-9.0+12.0*xi)*eta+(18.0-24.0*xi)*eta*eta+2.0*((6.0-8.0*xi)*eta+(-12.0+16.0*xi)*eta*eta)*zeta;
return_value[1][0] = 3.0-4.0*xi+2.0*(-6.0+8.0*xi)*eta+(-9.0+12.0*xi+2.0*(18.0-24.0*xi)*eta)*zeta+(6.0-8.0*xi+2.0*(-12.0+16.0*xi)*eta)*zeta*zeta;
return_value[1][1] = 4.0-12.0*xi+8.0*xi*xi+(-12.0+36.0*xi-24.0*xi*xi)*zeta+(8.0-24.0*xi+16.0*xi*xi)*zeta*zeta;
return_value[1][2] = 3.0-9.0*xi+6.0*xi*xi+2.0*(-6.0+18.0*xi-12.0*xi*xi)*eta+2.0*(-2.0+6.0*xi-4.0*xi*xi+2.0*(4.0-12.0*xi+8.0*xi*xi)*eta)*zeta;
return_value[2][0] = (-9.0+12.0*xi)*eta+(18.0-24.0*xi)*eta*eta+2.0*((6.0-8.0*xi)*eta+(-12.0+16.0*xi)*eta*eta)*zeta;
return_value[2][1] = 3.0-9.0*xi+6.0*xi*xi+2.0*(-6.0+18.0*xi-12.0*xi*xi)*eta+2.0*(-2.0+6.0*xi-4.0*xi*xi+2.0*(4.0-12.0*xi+8.0*xi*xi)*eta)*zeta;
return_value[2][2] = 2.0*(-2.0+6.0*xi-4.0*xi*xi)*eta+2.0*(4.0-12.0*xi+8.0*xi*xi)*eta*eta;
break;
case 5:
return_value[0][0] = -4.0*eta+8.0*eta*eta+(12.0*eta-24.0*eta*eta)*zeta+(-8.0*eta+16.0*eta*eta)*zeta*zeta;
return_value[0][1] = 1.0-4.0*xi+2.0*(-2.0+8.0*xi)*eta+(-3.0+12.0*xi+2.0*(6.0-24.0*xi)*eta)*zeta+(2.0-8.0*xi+2.0*(-4.0+16.0*xi)*eta)*zeta*zeta;
return_value[0][2] = (-3.0+12.0*xi)*eta+(6.0-24.0*xi)*eta*eta+2.0*((2.0-8.0*xi)*eta+(-4.0+16.0*xi)*eta*eta)*zeta;
return_value[1][0] = 1.0-4.0*xi+2.0*(-2.0+8.0*xi)*eta+(-3.0+12.0*xi+2.0*(6.0-24.0*xi)*eta)*zeta+(2.0-8.0*xi+2.0*(-4.0+16.0*xi)*eta)*zeta*zeta;
return_value[1][1] = -4.0*xi+8.0*xi*xi+(12.0*xi-24.0*xi*xi)*zeta+(-8.0*xi+16.0*xi*xi)*zeta*zeta;
return_value[1][2] = -3.0*xi+6.0*xi*xi+2.0*(6.0*xi-12.0*xi*xi)*eta+2.0*(2.0*xi-4.0*xi*xi+2.0*(-4.0*xi+8.0*xi*xi)*eta)*zeta;
return_value[2][0] = (-3.0+12.0*xi)*eta+(6.0-24.0*xi)*eta*eta+2.0*((2.0-8.0*xi)*eta+(-4.0+16.0*xi)*eta*eta)*zeta;
return_value[2][1] = -3.0*xi+6.0*xi*xi+2.0*(6.0*xi-12.0*xi*xi)*eta+2.0*(2.0*xi-4.0*xi*xi+2.0*(-4.0*xi+8.0*xi*xi)*eta)*zeta;
return_value[2][2] = 2.0*(2.0*xi-4.0*xi*xi)*eta+2.0*(-4.0*xi+8.0*xi*xi)*eta*eta;
break;
case 6:
return_value[0][0] = (4.0*eta-8.0*eta*eta)*zeta+(-8.0*eta+16.0*eta*eta)*zeta*zeta;
return_value[0][1] = (-1.0+4.0*xi+2.0*(2.0-8.0*xi)*eta)*zeta+(2.0-8.0*xi+2.0*(-4.0+16.0*xi)*eta)*zeta*zeta;
return_value[0][2] = (-1.0+4.0*xi)*eta+(2.0-8.0*xi)*eta*eta+2.0*((2.0-8.0*xi)*eta+(-4.0+16.0*xi)*eta*eta)*zeta;
return_value[1][0] = (-1.0+4.0*xi+2.0*(2.0-8.0*xi)*eta)*zeta+(2.0-8.0*xi+2.0*(-4.0+16.0*xi)*eta)*zeta*zeta;
return_value[1][1] = (4.0*xi-8.0*xi*xi)*zeta+(-8.0*xi+16.0*xi*xi)*zeta*zeta;
return_value[1][2] = -xi+2.0*xi*xi+2.0*(2.0*xi-4.0*xi*xi)*eta+2.0*(2.0*xi-4.0*xi*xi+2.0*(-4.0*xi+8.0*xi*xi)*eta)*zeta;
return_value[2][0] = (-1.0+4.0*xi)*eta+(2.0-8.0*xi)*eta*eta+2.0*((2.0-8.0*xi)*eta+(-4.0+16.0*xi)*eta*eta)*zeta;
return_value[2][1] = -xi+2.0*xi*xi+2.0*(2.0*xi-4.0*xi*xi)*eta+2.0*(2.0*xi-4.0*xi*xi+2.0*(-4.0*xi+8.0*xi*xi)*eta)*zeta;
return_value[2][2] = 2.0*(2.0*xi-4.0*xi*xi)*eta+2.0*(-4.0*xi+8.0*xi*xi)*eta*eta;
break;
case 7:
return_value[0][0] = (4.0*eta-8.0*eta*eta)*zeta+(-8.0*eta+16.0*eta*eta)*zeta*zeta;
return_value[0][1] = (-3.0+4.0*xi+2.0*(6.0-8.0*xi)*eta)*zeta+(6.0-8.0*xi+2.0*(-12.0+16.0*xi)*eta)*zeta*zeta;
return_value[0][2] = (-3.0+4.0*xi)*eta+(6.0-8.0*xi)*eta*eta+2.0*((6.0-8.0*xi)*eta+(-12.0+16.0*xi)*eta*eta)*zeta;
return_value[1][0] = (-3.0+4.0*xi+2.0*(6.0-8.0*xi)*eta)*zeta+(6.0-8.0*xi+2.0*(-12.0+16.0*xi)*eta)*zeta*zeta;
return_value[1][1] = (-4.0+12.0*xi-8.0*xi*xi)*zeta+(8.0-24.0*xi+16.0*xi*xi)*zeta*zeta;
return_value[1][2] = 1.0-3.0*xi+2.0*xi*xi+2.0*(-2.0+6.0*xi-4.0*xi*xi)*eta+2.0*(-2.0+6.0*xi-4.0*xi*xi+2.0*(4.0-12.0*xi+8.0*xi*xi)*eta)*zeta;
return_value[2][0] = (-3.0+4.0*xi)*eta+(6.0-8.0*xi)*eta*eta+2.0*((6.0-8.0*xi)*eta+(-12.0+16.0*xi)*eta*eta)*zeta;
return_value[2][1] = 1.0-3.0*xi+2.0*xi*xi+2.0*(-2.0+6.0*xi-4.0*xi*xi)*eta+2.0*(-2.0+6.0*xi-4.0*xi*xi+2.0*(4.0-12.0*xi+8.0*xi*xi)*eta)*zeta;
return_value[2][2] = 2.0*(-2.0+6.0*xi-4.0*xi*xi)*eta+2.0*(4.0-12.0*xi+8.0*xi*xi)*eta*eta;
break;
case 8:
return_value[0][0] = -8.0+24.0*eta-16.0*eta*eta+(24.0-72.0*eta+48.0*eta*eta)*zeta+(-16.0+48.0*eta-32.0*eta*eta)*zeta*zeta;
return_value[0][1] = -12.0+24.0*xi+2.0*(8.0-16.0*xi)*eta+(36.0-72.0*xi+2.0*(-24.0+48.0*xi)*eta)*zeta+(-24.0+48.0*xi+2.0*(16.0-32.0*xi)*eta)*zeta*zeta;
return_value[0][2] = -12.0+24.0*xi+(36.0-72.0*xi)*eta+(-24.0+48.0*xi)*eta*eta+2.0*(8.0-16.0*xi+(-24.0+48.0*xi)*eta+(16.0-32.0*xi)*eta*eta)*zeta;
return_value[1][0] = -12.0+24.0*xi+2.0*(8.0-16.0*xi)*eta+(36.0-72.0*xi+2.0*(-24.0+48.0*xi)*eta)*zeta+(-24.0+48.0*xi+2.0*(16.0-32.0*xi)*eta)*zeta*zeta;
return_value[1][1] = 16.0*xi-16.0*xi*xi+(-48.0*xi+48.0*xi*xi)*zeta+(32.0*xi-32.0*xi*xi)*zeta*zeta;
return_value[1][2] = 36.0*xi-36.0*xi*xi+2.0*(-24.0*xi+24.0*xi*xi)*eta+2.0*(-24.0*xi+24.0*xi*xi+2.0*(16.0*xi-16.0*xi*xi)*eta)*zeta;
return_value[2][0] = -12.0+24.0*xi+(36.0-72.0*xi)*eta+(-24.0+48.0*xi)*eta*eta+2.0*(8.0-16.0*xi+(-24.0+48.0*xi)*eta+(16.0-32.0*xi)*eta*eta)*zeta;
return_value[2][1] = 36.0*xi-36.0*xi*xi+2.0*(-24.0*xi+24.0*xi*xi)*eta+2.0*(-24.0*xi+24.0*xi*xi+2.0*(16.0*xi-16.0*xi*xi)*eta)*zeta;
return_value[2][2] = 16.0*xi-16.0*xi*xi+2.0*(-24.0*xi+24.0*xi*xi)*eta+2.0*(16.0*xi-16.0*xi*xi)*eta*eta;
break;
case 9:
return_value[0][0] = (16.0-48.0*eta+32.0*eta*eta)*zeta+(-16.0+48.0*eta-32.0*eta*eta)*zeta*zeta;
return_value[0][1] = (12.0-48.0*xi+2.0*(-8.0+32.0*xi)*eta)*zeta+(-12.0+48.0*xi+2.0*(8.0-32.0*xi)*eta)*zeta*zeta;
return_value[0][2] = -4.0+16.0*xi+(12.0-48.0*xi)*eta+(-8.0+32.0*xi)*eta*eta+2.0*(4.0-16.0*xi+(-12.0+48.0*xi)*eta+(8.0-32.0*xi)*eta*eta)*zeta;
return_value[1][0] = (12.0-48.0*xi+2.0*(-8.0+32.0*xi)*eta)*zeta+(-12.0+48.0*xi+2.0*(8.0-32.0*xi)*eta)*zeta*zeta;
return_value[1][1] = (-16.0*xi+32.0*xi*xi)*zeta+(16.0*xi-32.0*xi*xi)*zeta*zeta;
return_value[1][2] = 12.0*xi-24.0*xi*xi+2.0*(-8.0*xi+16.0*xi*xi)*eta+2.0*(-12.0*xi+24.0*xi*xi+2.0*(8.0*xi-16.0*xi*xi)*eta)*zeta;
return_value[2][0] = -4.0+16.0*xi+(12.0-48.0*xi)*eta+(-8.0+32.0*xi)*eta*eta+2.0*(4.0-16.0*xi+(-12.0+48.0*xi)*eta+(8.0-32.0*xi)*eta*eta)*zeta;
return_value[2][1] = 12.0*xi-24.0*xi*xi+2.0*(-8.0*xi+16.0*xi*xi)*eta+2.0*(-12.0*xi+24.0*xi*xi+2.0*(8.0*xi-16.0*xi*xi)*eta)*zeta;
return_value[2][2] = 8.0*xi-16.0*xi*xi+2.0*(-12.0*xi+24.0*xi*xi)*eta+2.0*(8.0*xi-16.0*xi*xi)*eta*eta;
return_value[0][0] = (8.0-24.0*eta+16.0*eta*eta)*zeta+(-16.0+48.0*eta-32.0*eta*eta)*zeta*zeta;
return_value[0][1] = (12.0-24.0*xi+2.0*(-8.0+16.0*xi)*eta)*zeta+(-24.0+48.0*xi+2.0*(16.0-32.0*xi)*eta)*zeta*zeta;
return_value[0][2] = -4.0+8.0*xi+(12.0-24.0*xi)*eta+(-8.0+16.0*xi)*eta*eta+2.0*(8.0-16.0*xi+(-24.0+48.0*xi)*eta+(16.0-32.0*xi)*eta*eta)*zeta;
return_value[1][0] = (12.0-24.0*xi+2.0*(-8.0+16.0*xi)*eta)*zeta+(-24.0+48.0*xi+2.0*(16.0-32.0*xi)*eta)*zeta*zeta;
return_value[1][1] = (-16.0*xi+16.0*xi*xi)*zeta+(32.0*xi-32.0*xi*xi)*zeta*zeta;
return_value[1][2] = 12.0*xi-12.0*xi*xi+2.0*(-8.0*xi+8.0*xi*xi)*eta+2.0*(-24.0*xi+24.0*xi*xi+2.0*(16.0*xi-16.0*xi*xi)*eta)*zeta;
return_value[2][0] = -4.0+8.0*xi+(12.0-24.0*xi)*eta+(-8.0+16.0*xi)*eta*eta+2.0*(8.0-16.0*xi+(-24.0+48.0*xi)*eta+(16.0-32.0*xi)*eta*eta)*zeta;
return_value[2][1] = 12.0*xi-12.0*xi*xi+2.0*(-8.0*xi+8.0*xi*xi)*eta+2.0*(-24.0*xi+24.0*xi*xi+2.0*(16.0*xi-16.0*xi*xi)*eta)*zeta;
return_value[2][2] = 16.0*xi-16.0*xi*xi+2.0*(-24.0*xi+24.0*xi*xi)*eta+2.0*(16.0*xi-16.0*xi*xi)*eta*eta;
return_value[0][0] = (16.0-48.0*eta+32.0*eta*eta)*zeta+(-16.0+48.0*eta-32.0*eta*eta)*zeta*zeta;
return_value[0][1] = (36.0-48.0*xi+2.0*(-24.0+32.0*xi)*eta)*zeta+(-36.0+48.0*xi+2.0*(24.0-32.0*xi)*eta)*zeta*zeta;
return_value[0][2] = -12.0+16.0*xi+(36.0-48.0*xi)*eta+(-24.0+32.0*xi)*eta*eta+2.0*(12.0-16.0*xi+(-36.0+48.0*xi)*eta+(24.0-32.0*xi)*eta*eta)*zeta;
return_value[1][0] = (36.0-48.0*xi+2.0*(-24.0+32.0*xi)*eta)*zeta+(-36.0+48.0*xi+2.0*(24.0-32.0*xi)*eta)*zeta*zeta;
return_value[1][1] = (16.0-48.0*xi+32.0*xi*xi)*zeta+(-16.0+48.0*xi-32.0*xi*xi)*zeta*zeta;
return_value[1][2] = -12.0+36.0*xi-24.0*xi*xi+2.0*(8.0-24.0*xi+16.0*xi*xi)*eta+2.0*(12.0-36.0*xi+24.0*xi*xi+2.0*(-8.0+24.0*xi-16.0*xi*xi)*eta)*zeta;
return_value[2][0] = -12.0+16.0*xi+(36.0-48.0*xi)*eta+(-24.0+32.0*xi)*eta*eta+2.0*(12.0-16.0*xi+(-36.0+48.0*xi)*eta+(24.0-32.0*xi)*eta*eta)*zeta;
return_value[2][1] = -12.0+36.0*xi-24.0*xi*xi+2.0*(8.0-24.0*xi+16.0*xi*xi)*eta+2.0*(12.0-36.0*xi+24.0*xi*xi+2.0*(-8.0+24.0*xi-16.0*xi*xi)*eta)*zeta;
return_value[2][2] = -8.0+24.0*xi-16.0*xi*xi+2.0*(12.0-36.0*xi+24.0*xi*xi)*eta+2.0*(-8.0+24.0*xi-16.0*xi*xi)*eta*eta;
return_value[0][0] = 8.0*eta-16.0*eta*eta+(-24.0*eta+48.0*eta*eta)*zeta+(16.0*eta-32.0*eta*eta)*zeta*zeta;
return_value[0][1] = -4.0+8.0*xi+2.0*(8.0-16.0*xi)*eta+(12.0-24.0*xi+2.0*(-24.0+48.0*xi)*eta)*zeta+(-8.0+16.0*xi+2.0*(16.0-32.0*xi)*eta)*zeta*zeta;
return_value[0][2] = (12.0-24.0*xi)*eta+(-24.0+48.0*xi)*eta*eta+2.0*((-8.0+16.0*xi)*eta+(16.0-32.0*xi)*eta*eta)*zeta;
return_value[1][0] = -4.0+8.0*xi+2.0*(8.0-16.0*xi)*eta+(12.0-24.0*xi+2.0*(-24.0+48.0*xi)*eta)*zeta+(-8.0+16.0*xi+2.0*(16.0-32.0*xi)*eta)*zeta*zeta;
return_value[1][1] = 16.0*xi-16.0*xi*xi+(-48.0*xi+48.0*xi*xi)*zeta+(32.0*xi-32.0*xi*xi)*zeta*zeta;
return_value[1][2] = 12.0*xi-12.0*xi*xi+2.0*(-24.0*xi+24.0*xi*xi)*eta+2.0*(-8.0*xi+8.0*xi*xi+2.0*(16.0*xi-16.0*xi*xi)*eta)*zeta;
return_value[2][0] = (12.0-24.0*xi)*eta+(-24.0+48.0*xi)*eta*eta+2.0*((-8.0+16.0*xi)*eta+(16.0-32.0*xi)*eta*eta)*zeta;
return_value[2][1] = 12.0*xi-12.0*xi*xi+2.0*(-24.0*xi+24.0*xi*xi)*eta+2.0*(-8.0*xi+8.0*xi*xi+2.0*(16.0*xi-16.0*xi*xi)*eta)*zeta;
return_value[2][2] = 2.0*(-8.0*xi+8.0*xi*xi)*eta+2.0*(16.0*xi-16.0*xi*xi)*eta*eta;
return_value[0][0] = (-16.0*eta+32.0*eta*eta)*zeta+(16.0*eta-32.0*eta*eta)*zeta*zeta;
return_value[0][1] = (4.0-16.0*xi+2.0*(-8.0+32.0*xi)*eta)*zeta+(-4.0+16.0*xi+2.0*(8.0-32.0*xi)*eta)*zeta*zeta;
return_value[0][2] = (4.0-16.0*xi)*eta+(-8.0+32.0*xi)*eta*eta+2.0*((-4.0+16.0*xi)*eta+(8.0-32.0*xi)*eta*eta)*zeta;
return_value[1][0] = (4.0-16.0*xi+2.0*(-8.0+32.0*xi)*eta)*zeta+(-4.0+16.0*xi+2.0*(8.0-32.0*xi)*eta)*zeta*zeta;
return_value[1][1] = (-16.0*xi+32.0*xi*xi)*zeta+(16.0*xi-32.0*xi*xi)*zeta*zeta;
return_value[1][2] = 4.0*xi-8.0*xi*xi+2.0*(-8.0*xi+16.0*xi*xi)*eta+2.0*(-4.0*xi+8.0*xi*xi+2.0*(8.0*xi-16.0*xi*xi)*eta)*zeta;
return_value[2][0] = (4.0-16.0*xi)*eta+(-8.0+32.0*xi)*eta*eta+2.0*((-4.0+16.0*xi)*eta+(8.0-32.0*xi)*eta*eta)*zeta;
return_value[2][1] = 4.0*xi-8.0*xi*xi+2.0*(-8.0*xi+16.0*xi*xi)*eta+2.0*(-4.0*xi+8.0*xi*xi+2.0*(8.0*xi-16.0*xi*xi)*eta)*zeta;
return_value[2][2] = 2.0*(-4.0*xi+8.0*xi*xi)*eta+2.0*(8.0*xi-16.0*xi*xi)*eta*eta;
return_value[0][0] = (-8.0*eta+16.0*eta*eta)*zeta+(16.0*eta-32.0*eta*eta)*zeta*zeta;
return_value[0][1] = (4.0-8.0*xi+2.0*(-8.0+16.0*xi)*eta)*zeta+(-8.0+16.0*xi+2.0*(16.0-32.0*xi)*eta)*zeta*zeta;
return_value[0][2] = (4.0-8.0*xi)*eta+(-8.0+16.0*xi)*eta*eta+2.0*((-8.0+16.0*xi)*eta+(16.0-32.0*xi)*eta*eta)*zeta;
return_value[1][0] = (4.0-8.0*xi+2.0*(-8.0+16.0*xi)*eta)*zeta+(-8.0+16.0*xi+2.0*(16.0-32.0*xi)*eta)*zeta*zeta;
return_value[1][1] = (-16.0*xi+16.0*xi*xi)*zeta+(32.0*xi-32.0*xi*xi)*zeta*zeta;
return_value[1][2] = 4.0*xi-4.0*xi*xi+2.0*(-8.0*xi+8.0*xi*xi)*eta+2.0*(-8.0*xi+8.0*xi*xi+2.0*(16.0*xi-16.0*xi*xi)*eta)*zeta;
return_value[2][0] = (4.0-8.0*xi)*eta+(-8.0+16.0*xi)*eta*eta+2.0*((-8.0+16.0*xi)*eta+(16.0-32.0*xi)*eta*eta)*zeta;
return_value[2][1] = 4.0*xi-4.0*xi*xi+2.0*(-8.0*xi+8.0*xi*xi)*eta+2.0*(-8.0*xi+8.0*xi*xi+2.0*(16.0*xi-16.0*xi*xi)*eta)*zeta;
return_value[2][2] = 2.0*(-8.0*xi+8.0*xi*xi)*eta+2.0*(16.0*xi-16.0*xi*xi)*eta*eta;
return_value[0][0] = (-16.0*eta+32.0*eta*eta)*zeta+(16.0*eta-32.0*eta*eta)*zeta*zeta;
return_value[0][1] = (12.0-16.0*xi+2.0*(-24.0+32.0*xi)*eta)*zeta+(-12.0+16.0*xi+2.0*(24.0-32.0*xi)*eta)*zeta*zeta;
return_value[0][2] = (12.0-16.0*xi)*eta+(-24.0+32.0*xi)*eta*eta+2.0*((-12.0+16.0*xi)*eta+(24.0-32.0*xi)*eta*eta)*zeta;
return_value[1][0] = (12.0-16.0*xi+2.0*(-24.0+32.0*xi)*eta)*zeta+(-12.0+16.0*xi+2.0*(24.0-32.0*xi)*eta)*zeta*zeta;
return_value[1][1] = (16.0-48.0*xi+32.0*xi*xi)*zeta+(-16.0+48.0*xi-32.0*xi*xi)*zeta*zeta;
return_value[1][2] = -4.0+12.0*xi-8.0*xi*xi+2.0*(8.0-24.0*xi+16.0*xi*xi)*eta+2.0*(4.0-12.0*xi+8.0*xi*xi+2.0*(-8.0+24.0*xi-16.0*xi*xi)*eta)*zeta;
return_value[2][0] = (12.0-16.0*xi)*eta+(-24.0+32.0*xi)*eta*eta+2.0*((-12.0+16.0*xi)*eta+(24.0-32.0*xi)*eta*eta)*zeta;
return_value[2][1] = -4.0+12.0*xi-8.0*xi*xi+2.0*(8.0-24.0*xi+16.0*xi*xi)*eta+2.0*(4.0-12.0*xi+8.0*xi*xi+2.0*(-8.0+24.0*xi-16.0*xi*xi)*eta)*zeta;
return_value[2][2] = 2.0*(4.0-12.0*xi+8.0*xi*xi)*eta+2.0*(-8.0+24.0*xi-16.0*xi*xi)*eta*eta;
return_value[0][0] = 16.0*eta-16.0*eta*eta+(-48.0*eta+48.0*eta*eta)*zeta+(32.0*eta-32.0*eta*eta)*zeta*zeta;
return_value[0][1] = -12.0+16.0*xi+2.0*(12.0-16.0*xi)*eta+(36.0-48.0*xi+2.0*(-36.0+48.0*xi)*eta)*zeta+(-24.0+32.0*xi+2.0*(24.0-32.0*xi)*eta)*zeta*zeta;
return_value[0][2] = (36.0-48.0*xi)*eta+(-36.0+48.0*xi)*eta*eta+2.0*((-24.0+32.0*xi)*eta+(24.0-32.0*xi)*eta*eta)*zeta;
return_value[1][0] = -12.0+16.0*xi+2.0*(12.0-16.0*xi)*eta+(36.0-48.0*xi+2.0*(-36.0+48.0*xi)*eta)*zeta+(-24.0+32.0*xi+2.0*(24.0-32.0*xi)*eta)*zeta*zeta;
return_value[1][1] = -8.0+24.0*xi-16.0*xi*xi+(24.0-72.0*xi+48.0*xi*xi)*zeta+(-16.0+48.0*xi-32.0*xi*xi)*zeta*zeta;
return_value[1][2] = -12.0+36.0*xi-24.0*xi*xi+2.0*(12.0-36.0*xi+24.0*xi*xi)*eta+2.0*(8.0-24.0*xi+16.0*xi*xi+2.0*(-8.0+24.0*xi-16.0*xi*xi)*eta)*zeta;
return_value[2][0] = (36.0-48.0*xi)*eta+(-36.0+48.0*xi)*eta*eta+2.0*((-24.0+32.0*xi)*eta+(24.0-32.0*xi)*eta*eta)*zeta;
return_value[2][1] = -12.0+36.0*xi-24.0*xi*xi+2.0*(12.0-36.0*xi+24.0*xi*xi)*eta+2.0*(8.0-24.0*xi+16.0*xi*xi+2.0*(-8.0+24.0*xi-16.0*xi*xi)*eta)*zeta;
return_value[2][2] = 2.0*(8.0-24.0*xi+16.0*xi*xi)*eta+2.0*(-8.0+24.0*xi-16.0*xi*xi)*eta*eta;
return_value[0][0] = 16.0*eta-16.0*eta*eta+(-48.0*eta+48.0*eta*eta)*zeta+(32.0*eta-32.0*eta*eta)*zeta*zeta;
return_value[0][1] = -4.0+16.0*xi+2.0*(4.0-16.0*xi)*eta+(12.0-48.0*xi+2.0*(-12.0+48.0*xi)*eta)*zeta+(-8.0+32.0*xi+2.0*(8.0-32.0*xi)*eta)*zeta*zeta;
return_value[0][2] = (12.0-48.0*xi)*eta+(-12.0+48.0*xi)*eta*eta+2.0*((-8.0+32.0*xi)*eta+(8.0-32.0*xi)*eta*eta)*zeta;
return_value[1][0] = -4.0+16.0*xi+2.0*(4.0-16.0*xi)*eta+(12.0-48.0*xi+2.0*(-12.0+48.0*xi)*eta)*zeta+(-8.0+32.0*xi+2.0*(8.0-32.0*xi)*eta)*zeta*zeta;
return_value[1][1] = 8.0*xi-16.0*xi*xi+(-24.0*xi+48.0*xi*xi)*zeta+(16.0*xi-32.0*xi*xi)*zeta*zeta;
return_value[1][2] = 12.0*xi-24.0*xi*xi+2.0*(-12.0*xi+24.0*xi*xi)*eta+2.0*(-8.0*xi+16.0*xi*xi+2.0*(8.0*xi-16.0*xi*xi)*eta)*zeta;
return_value[2][0] = (12.0-48.0*xi)*eta+(-12.0+48.0*xi)*eta*eta+2.0*((-8.0+32.0*xi)*eta+(8.0-32.0*xi)*eta*eta)*zeta;
return_value[2][1] = 12.0*xi-24.0*xi*xi+2.0*(-12.0*xi+24.0*xi*xi)*eta+2.0*(-8.0*xi+16.0*xi*xi+2.0*(8.0*xi-16.0*xi*xi)*eta)*zeta;
return_value[2][2] = 2.0*(-8.0*xi+16.0*xi*xi)*eta+2.0*(8.0*xi-16.0*xi*xi)*eta*eta;
return_value[0][0] = (-16.0*eta+16.0*eta*eta)*zeta+(32.0*eta-32.0*eta*eta)*zeta*zeta;
return_value[0][1] = (4.0-16.0*xi+2.0*(-4.0+16.0*xi)*eta)*zeta+(-8.0+32.0*xi+2.0*(8.0-32.0*xi)*eta)*zeta*zeta;
return_value[0][2] = (4.0-16.0*xi)*eta+(-4.0+16.0*xi)*eta*eta+2.0*((-8.0+32.0*xi)*eta+(8.0-32.0*xi)*eta*eta)*zeta;
return_value[1][0] = (4.0-16.0*xi+2.0*(-4.0+16.0*xi)*eta)*zeta+(-8.0+32.0*xi+2.0*(8.0-32.0*xi)*eta)*zeta*zeta;
return_value[1][1] = (-8.0*xi+16.0*xi*xi)*zeta+(16.0*xi-32.0*xi*xi)*zeta*zeta;
return_value[1][2] = 4.0*xi-8.0*xi*xi+2.0*(-4.0*xi+8.0*xi*xi)*eta+2.0*(-8.0*xi+16.0*xi*xi+2.0*(8.0*xi-16.0*xi*xi)*eta)*zeta;
return_value[2][0] = (4.0-16.0*xi)*eta+(-4.0+16.0*xi)*eta*eta+2.0*((-8.0+32.0*xi)*eta+(8.0-32.0*xi)*eta*eta)*zeta;
return_value[2][1] = 4.0*xi-8.0*xi*xi+2.0*(-4.0*xi+8.0*xi*xi)*eta+2.0*(-8.0*xi+16.0*xi*xi+2.0*(8.0*xi-16.0*xi*xi)*eta)*zeta;
return_value[2][2] = 2.0*(-8.0*xi+16.0*xi*xi)*eta+2.0*(8.0*xi-16.0*xi*xi)*eta*eta;
return_value[0][0] = (-16.0*eta+16.0*eta*eta)*zeta+(32.0*eta-32.0*eta*eta)*zeta*zeta;
return_value[0][1] = (12.0-16.0*xi+2.0*(-12.0+16.0*xi)*eta)*zeta+(-24.0+32.0*xi+2.0*(24.0-32.0*xi)*eta)*zeta*zeta;
return_value[0][2] = (12.0-16.0*xi)*eta+(-12.0+16.0*xi)*eta*eta+2.0*((-24.0+32.0*xi)*eta+(24.0-32.0*xi)*eta*eta)*zeta;
return_value[1][0] = (12.0-16.0*xi+2.0*(-12.0+16.0*xi)*eta)*zeta+(-24.0+32.0*xi+2.0*(24.0-32.0*xi)*eta)*zeta*zeta;
return_value[1][1] = (8.0-24.0*xi+16.0*xi*xi)*zeta+(-16.0+48.0*xi-32.0*xi*xi)*zeta*zeta;
return_value[1][2] = -4.0+12.0*xi-8.0*xi*xi+2.0*(4.0-12.0*xi+8.0*xi*xi)*eta+2.0*(8.0-24.0*xi+16.0*xi*xi+2.0*(-8.0+24.0*xi-16.0*xi*xi)*eta)*zeta;
return_value[2][0] = (12.0-16.0*xi)*eta+(-12.0+16.0*xi)*eta*eta+2.0*((-24.0+32.0*xi)*eta+(24.0-32.0*xi)*eta*eta)*zeta;
return_value[2][1] = -4.0+12.0*xi-8.0*xi*xi+2.0*(4.0-12.0*xi+8.0*xi*xi)*eta+2.0*(8.0-24.0*xi+16.0*xi*xi+2.0*(-8.0+24.0*xi-16.0*xi*xi)*eta)*zeta;
return_value[2][2] = 2.0*(8.0-24.0*xi+16.0*xi*xi)*eta+2.0*(-8.0+24.0*xi-16.0*xi*xi)*eta*eta;
return_value[0][0] = (-32.0+96.0*eta-64.0*eta*eta)*zeta+(32.0-96.0*eta+64.0*eta*eta)*zeta*zeta;
return_value[0][1] = (-48.0+96.0*xi+2.0*(32.0-64.0*xi)*eta)*zeta+(48.0-96.0*xi+2.0*(-32.0+64.0*xi)*eta)*zeta*zeta;
return_value[0][2] = 16.0-32.0*xi+(-48.0+96.0*xi)*eta+(32.0-64.0*xi)*eta*eta+2.0*(-16.0+32.0*xi+(48.0-96.0*xi)*eta+(-32.0+64.0*xi)*eta*eta)*zeta;
return_value[1][0] = (-48.0+96.0*xi+2.0*(32.0-64.0*xi)*eta)*zeta+(48.0-96.0*xi+2.0*(-32.0+64.0*xi)*eta)*zeta*zeta;
return_value[1][1] = (64.0*xi-64.0*xi*xi)*zeta+(-64.0*xi+64.0*xi*xi)*zeta*zeta;
return_value[1][2] = -48.0*xi+48.0*xi*xi+2.0*(32.0*xi-32.0*xi*xi)*eta+2.0*(48.0*xi-48.0*xi*xi+2.0*(-32.0*xi+32.0*xi*xi)*eta)*zeta;
return_value[2][0] = 16.0-32.0*xi+(-48.0+96.0*xi)*eta+(32.0-64.0*xi)*eta*eta+2.0*(-16.0+32.0*xi+(48.0-96.0*xi)*eta+(-32.0+64.0*xi)*eta*eta)*zeta;
return_value[2][1] = -48.0*xi+48.0*xi*xi+2.0*(32.0*xi-32.0*xi*xi)*eta+2.0*(48.0*xi-48.0*xi*xi+2.0*(-32.0*xi+32.0*xi*xi)*eta)*zeta;
return_value[2][2] = -32.0*xi+32.0*xi*xi+2.0*(48.0*xi-48.0*xi*xi)*eta+2.0*(-32.0*xi+32.0*xi*xi)*eta*eta;
return_value[0][0] = (32.0*eta-64.0*eta*eta)*zeta+(-32.0*eta+64.0*eta*eta)*zeta*zeta;
return_value[0][1] = (-16.0+32.0*xi+2.0*(32.0-64.0*xi)*eta)*zeta+(16.0-32.0*xi+2.0*(-32.0+64.0*xi)*eta)*zeta*zeta;
return_value[0][2] = (-16.0+32.0*xi)*eta+(32.0-64.0*xi)*eta*eta+2.0*((16.0-32.0*xi)*eta+(-32.0+64.0*xi)*eta*eta)*zeta;
return_value[1][0] = (-16.0+32.0*xi+2.0*(32.0-64.0*xi)*eta)*zeta+(16.0-32.0*xi+2.0*(-32.0+64.0*xi)*eta)*zeta*zeta;
return_value[1][1] = (64.0*xi-64.0*xi*xi)*zeta+(-64.0*xi+64.0*xi*xi)*zeta*zeta;
return_value[1][2] = -16.0*xi+16.0*xi*xi+2.0*(32.0*xi-32.0*xi*xi)*eta+2.0*(16.0*xi-16.0*xi*xi+2.0*(-32.0*xi+32.0*xi*xi)*eta)*zeta;
return_value[2][0] = (-16.0+32.0*xi)*eta+(32.0-64.0*xi)*eta*eta+2.0*((16.0-32.0*xi)*eta+(-32.0+64.0*xi)*eta*eta)*zeta;
return_value[2][1] = -16.0*xi+16.0*xi*xi+2.0*(32.0*xi-32.0*xi*xi)*eta+2.0*(16.0*xi-16.0*xi*xi+2.0*(-32.0*xi+32.0*xi*xi)*eta)*zeta;
return_value[2][2] = 2.0*(16.0*xi-16.0*xi*xi)*eta+2.0*(-32.0*xi+32.0*xi*xi)*eta*eta;
return_value[0][0] = -32.0*eta+32.0*eta*eta+(96.0*eta-96.0*eta*eta)*zeta+(-64.0*eta+64.0*eta*eta)*zeta*zeta;
return_value[0][1] = 16.0-32.0*xi+2.0*(-16.0+32.0*xi)*eta+(-48.0+96.0*xi+2.0*(48.0-96.0*xi)*eta)*zeta+(32.0-64.0*xi+2.0*(-32.0+64.0*xi)*eta)*zeta*zeta;
return_value[0][2] = (-48.0+96.0*xi)*eta+(48.0-96.0*xi)*eta*eta+2.0*((32.0-64.0*xi)*eta+(-32.0+64.0*xi)*eta*eta)*zeta;
return_value[1][0] = 16.0-32.0*xi+2.0*(-16.0+32.0*xi)*eta+(-48.0+96.0*xi+2.0*(48.0-96.0*xi)*eta)*zeta+(32.0-64.0*xi+2.0*(-32.0+64.0*xi)*eta)*zeta*zeta;
return_value[1][1] = -32.0*xi+32.0*xi*xi+(96.0*xi-96.0*xi*xi)*zeta+(-64.0*xi+64.0*xi*xi)*zeta*zeta;
return_value[1][2] = -48.0*xi+48.0*xi*xi+2.0*(48.0*xi-48.0*xi*xi)*eta+2.0*(32.0*xi-32.0*xi*xi+2.0*(-32.0*xi+32.0*xi*xi)*eta)*zeta;
return_value[2][0] = (-48.0+96.0*xi)*eta+(48.0-96.0*xi)*eta*eta+2.0*((32.0-64.0*xi)*eta+(-32.0+64.0*xi)*eta*eta)*zeta;
return_value[2][1] = -48.0*xi+48.0*xi*xi+2.0*(48.0*xi-48.0*xi*xi)*eta+2.0*(32.0*xi-32.0*xi*xi+2.0*(-32.0*xi+32.0*xi*xi)*eta)*zeta;
return_value[2][2] = 2.0*(32.0*xi-32.0*xi*xi)*eta+2.0*(-32.0*xi+32.0*xi*xi)*eta*eta;
return_value[0][0] = (64.0*eta-64.0*eta*eta)*zeta+(-64.0*eta+64.0*eta*eta)*zeta*zeta;
return_value[0][1] = (-16.0+64.0*xi+2.0*(16.0-64.0*xi)*eta)*zeta+(16.0-64.0*xi+2.0*(-16.0+64.0*xi)*eta)*zeta*zeta;
return_value[0][2] = (-16.0+64.0*xi)*eta+(16.0-64.0*xi)*eta*eta+2.0*((16.0-64.0*xi)*eta+(-16.0+64.0*xi)*eta*eta)*zeta;
return_value[1][0] = (-16.0+64.0*xi+2.0*(16.0-64.0*xi)*eta)*zeta+(16.0-64.0*xi+2.0*(-16.0+64.0*xi)*eta)*zeta*zeta;
return_value[1][1] = (32.0*xi-64.0*xi*xi)*zeta+(-32.0*xi+64.0*xi*xi)*zeta*zeta;
return_value[1][2] = -16.0*xi+32.0*xi*xi+2.0*(16.0*xi-32.0*xi*xi)*eta+2.0*(16.0*xi-32.0*xi*xi+2.0*(-16.0*xi+32.0*xi*xi)*eta)*zeta;
return_value[2][0] = (-16.0+64.0*xi)*eta+(16.0-64.0*xi)*eta*eta+2.0*((16.0-64.0*xi)*eta+(-16.0+64.0*xi)*eta*eta)*zeta;
return_value[2][1] = -16.0*xi+32.0*xi*xi+2.0*(16.0*xi-32.0*xi*xi)*eta+2.0*(16.0*xi-32.0*xi*xi+2.0*(-16.0*xi+32.0*xi*xi)*eta)*zeta;
return_value[2][2] = 2.0*(16.0*xi-32.0*xi*xi)*eta+2.0*(-16.0*xi+32.0*xi*xi)*eta*eta;
return_value[0][0] = (32.0*eta-32.0*eta*eta)*zeta+(-64.0*eta+64.0*eta*eta)*zeta*zeta;
return_value[0][1] = (-16.0+32.0*xi+2.0*(16.0-32.0*xi)*eta)*zeta+(32.0-64.0*xi+2.0*(-32.0+64.0*xi)*eta)*zeta*zeta;
return_value[0][2] = (-16.0+32.0*xi)*eta+(16.0-32.0*xi)*eta*eta+2.0*((32.0-64.0*xi)*eta+(-32.0+64.0*xi)*eta*eta)*zeta;
return_value[1][0] = (-16.0+32.0*xi+2.0*(16.0-32.0*xi)*eta)*zeta+(32.0-64.0*xi+2.0*(-32.0+64.0*xi)*eta)*zeta*zeta;
return_value[1][1] = (32.0*xi-32.0*xi*xi)*zeta+(-64.0*xi+64.0*xi*xi)*zeta*zeta;
return_value[1][2] = -16.0*xi+16.0*xi*xi+2.0*(16.0*xi-16.0*xi*xi)*eta+2.0*(32.0*xi-32.0*xi*xi+2.0*(-32.0*xi+32.0*xi*xi)*eta)*zeta;
return_value[2][0] = (-16.0+32.0*xi)*eta+(16.0-32.0*xi)*eta*eta+2.0*((32.0-64.0*xi)*eta+(-32.0+64.0*xi)*eta*eta)*zeta;
return_value[2][1] = -16.0*xi+16.0*xi*xi+2.0*(16.0*xi-16.0*xi*xi)*eta+2.0*(32.0*xi-32.0*xi*xi+2.0*(-32.0*xi+32.0*xi*xi)*eta)*zeta;
return_value[2][2] = 2.0*(32.0*xi-32.0*xi*xi)*eta+2.0*(-32.0*xi+32.0*xi*xi)*eta*eta;
return_value[0][0] = (64.0*eta-64.0*eta*eta)*zeta+(-64.0*eta+64.0*eta*eta)*zeta*zeta;
return_value[0][1] = (-48.0+64.0*xi+2.0*(48.0-64.0*xi)*eta)*zeta+(48.0-64.0*xi+2.0*(-48.0+64.0*xi)*eta)*zeta*zeta;
return_value[0][2] = (-48.0+64.0*xi)*eta+(48.0-64.0*xi)*eta*eta+2.0*((48.0-64.0*xi)*eta+(-48.0+64.0*xi)*eta*eta)*zeta;
return_value[1][0] = (-48.0+64.0*xi+2.0*(48.0-64.0*xi)*eta)*zeta+(48.0-64.0*xi+2.0*(-48.0+64.0*xi)*eta)*zeta*zeta;
return_value[1][1] = (-32.0+96.0*xi-64.0*xi*xi)*zeta+(32.0-96.0*xi+64.0*xi*xi)*zeta*zeta;
return_value[1][2] = 16.0-48.0*xi+32.0*xi*xi+2.0*(-16.0+48.0*xi-32.0*xi*xi)*eta+2.0*(-16.0+48.0*xi-32.0*xi*xi+2.0*(16.0-48.0*xi+32.0*xi*xi)*eta)*zeta;
return_value[2][0] = (-48.0+64.0*xi)*eta+(48.0-64.0*xi)*eta*eta+2.0*((48.0-64.0*xi)*eta+(-48.0+64.0*xi)*eta*eta)*zeta;
return_value[2][1] = 16.0-48.0*xi+32.0*xi*xi+2.0*(-16.0+48.0*xi-32.0*xi*xi)*eta+2.0*(-16.0+48.0*xi-32.0*xi*xi+2.0*(16.0-48.0*xi+32.0*xi*xi)*eta)*zeta;
return_value[2][2] = 2.0*(-16.0+48.0*xi-32.0*xi*xi)*eta+2.0*(16.0-48.0*xi+32.0*xi*xi)*eta*eta;
return_value[0][0] = (-128.0*eta+128.0*eta*eta)*zeta+(128.0*eta-128.0*eta*eta)*zeta*zeta;
return_value[0][1] = (64.0-128.0*xi+2.0*(-64.0+128.0*xi)*eta)*zeta+(-64.0+128.0*xi+2.0*(64.0-128.0*xi)*eta)*zeta*zeta;
return_value[0][2] = (64.0-128.0*xi)*eta+(-64.0+128.0*xi)*eta*eta+2.0*((-64.0+128.0*xi)*eta+(64.0-128.0*xi)*eta*eta)*zeta;
return_value[1][0] = (64.0-128.0*xi+2.0*(-64.0+128.0*xi)*eta)*zeta+(-64.0+128.0*xi+2.0*(64.0-128.0*xi)*eta)*zeta*zeta;
return_value[1][1] = (-128.0*xi+128.0*xi*xi)*zeta+(128.0*xi-128.0*xi*xi)*zeta*zeta;
return_value[1][2] = 64.0*xi-64.0*xi*xi+2.0*(-64.0*xi+64.0*xi*xi)*eta+2.0*(-64.0*xi+64.0*xi*xi+2.0*(64.0*xi-64.0*xi*xi)*eta)*zeta;
return_value[2][0] = (64.0-128.0*xi)*eta+(-64.0+128.0*xi)*eta*eta+2.0*((-64.0+128.0*xi)*eta+(64.0-128.0*xi)*eta*eta)*zeta;
return_value[2][1] = 64.0*xi-64.0*xi*xi+2.0*(-64.0*xi+64.0*xi*xi)*eta+2.0*(-64.0*xi+64.0*xi*xi+2.0*(64.0*xi-64.0*xi*xi)*eta)*zeta;
return_value[2][2] = 2.0*(-64.0*xi+64.0*xi*xi)*eta+2.0*(64.0*xi-64.0*xi*xi)*eta*eta;
 break;
 
    };
  return return_value;
};



template <>
void FEQuadraticSub<3>::get_local_mass_matrix (const DoFHandler<3>::cell_iterator &,
					       const Boundary<3> &,
					       dFMatrix &local_mass_matrix) const {
  Assert (local_mass_matrix.n() == total_dofs,
	  ExcWrongFieldDimension(local_mass_matrix.n(),total_dofs));
  Assert (local_mass_matrix.m() == total_dofs,
	  ExcWrongFieldDimension(local_mass_matrix.m(),total_dofs));

  AssertThrow (false, ExcComputationNotUseful(3));
};



template <>
void FEQuadraticSub<3>::get_unit_support_points (vector<Point<3> > &unit_points) const {
  Assert (unit_points.size() == total_dofs,
	  ExcWrongFieldDimension (unit_points.size(), total_dofs));
  
  unit_points[0] = Point<3>(0, 0, 0);
  unit_points[1] = Point<3>(1, 0, 0);
  unit_points[2] = Point<3>(1, 0, 1);
  unit_points[3] = Point<3>(0, 0, 1);
  unit_points[4] = Point<3>(0, 1, 0);
  unit_points[5] = Point<3>(1, 1, 0);
  unit_points[6] = Point<3>(1, 1, 1);
  unit_points[7] = Point<3>(0, 1, 1);
  unit_points[8] = Point<3>(1/2, 0, 0);
  unit_points[9] = Point<3>(1, 0, 1/2);
  unit_points[10] = Point<3>(1/2, 0, 1);
  unit_points[11] = Point<3>(0, 0, 1/2);
  unit_points[12] = Point<3>(1/2, 1, 0);
  unit_points[13] = Point<3>(1, 1, 1/2);
  unit_points[14] = Point<3>(1/2, 1, 1);
  unit_points[15] = Point<3>(0, 1, 1/2);
  unit_points[16] = Point<3>(0, 1/2, 0);
  unit_points[17] = Point<3>(1, 1/2, 0);
  unit_points[18] = Point<3>(1, 1/2, 1);
  unit_points[19] = Point<3>(0, 1/2, 1);
  unit_points[20] = Point<3>(1/2, 0, 1/2);
  unit_points[21] = Point<3>(1/2, 1, 1/2);
  unit_points[22] = Point<3>(1/2, 1/2, 0);
  unit_points[23] = Point<3>(1, 1/2, 1/2);
  unit_points[24] = Point<3>(1/2, 1/2, 1);
  unit_points[25] = Point<3>(0, 1/2, 1/2);
  unit_points[26] = Point<3>(1/2, 1/2, 1/2);
};

  

template <>
void FEQuadraticSub<3>::get_support_points (const typename DoFHandler<3>::cell_iterator &cell,
					    const Boundary<3>&,
					    vector<Point<3> >  &support_points) const {
  Assert (support_points.size() == total_dofs,
	  ExcWrongFieldDimension (support_points.size(), total_dofs));

  const Point<3> vertices[8] = { cell->vertex(0),
				 cell->vertex(1),
				 cell->vertex(2),
				 cell->vertex(3),
				 cell->vertex(4),
				 cell->vertex(5),
				 cell->vertex(6),
				 cell->vertex(7)  };
  
  support_points[0](0) = vertices[0](0);
  support_points[0](1) = vertices[0](1);
  support_points[0](2) = vertices[0](2);
  support_points[1](0) = vertices[1](0);
  support_points[1](1) = vertices[1](1);
  support_points[1](2) = vertices[1](2);
  support_points[2](0) = vertices[2](0);
  support_points[2](1) = vertices[2](1);
  support_points[2](2) = vertices[2](2);
  support_points[3](0) = vertices[3](0);
  support_points[3](1) = vertices[3](1);
  support_points[3](2) = vertices[3](2);
  support_points[4](0) = vertices[4](0);
  support_points[4](1) = vertices[4](1);
  support_points[4](2) = vertices[4](2);
  support_points[5](0) = vertices[5](0);
  support_points[5](1) = vertices[5](1);
  support_points[5](2) = vertices[5](2);
  support_points[6](0) = vertices[6](0);
  support_points[6](1) = vertices[6](1);
  support_points[6](2) = vertices[6](2);
  support_points[7](0) = vertices[7](0);
  support_points[7](1) = vertices[7](1);
  support_points[7](2) = vertices[7](2);
  support_points[8](0) = vertices[0](0)/2.0+vertices[1](0)/2.0;
  support_points[8](1) = vertices[0](1)/2.0+vertices[1](1)/2.0;
  support_points[8](2) = vertices[0](2)/2.0+vertices[1](2)/2.0;
  support_points[9](0) = vertices[1](0)/2.0+vertices[2](0)/2.0;
  support_points[9](1) = vertices[1](1)/2.0+vertices[2](1)/2.0;
  support_points[9](2) = vertices[1](2)/2.0+vertices[2](2)/2.0;
  support_points[10](0) = vertices[2](0)/2.0+vertices[3](0)/2.0;
  support_points[10](1) = vertices[2](1)/2.0+vertices[3](1)/2.0;
  support_points[10](2) = vertices[2](2)/2.0+vertices[3](2)/2.0;
  support_points[11](0) = vertices[0](0)/2.0+vertices[3](0)/2.0;
  support_points[11](1) = vertices[0](1)/2.0+vertices[3](1)/2.0;
  support_points[11](2) = vertices[0](2)/2.0+vertices[3](2)/2.0;
  support_points[12](0) = vertices[4](0)/2.0+vertices[5](0)/2.0;
  support_points[12](1) = vertices[4](1)/2.0+vertices[5](1)/2.0;
  support_points[12](2) = vertices[4](2)/2.0+vertices[5](2)/2.0;
  support_points[13](0) = vertices[5](0)/2.0+vertices[6](0)/2.0;
  support_points[13](1) = vertices[5](1)/2.0+vertices[6](1)/2.0;
  support_points[13](2) = vertices[5](2)/2.0+vertices[6](2)/2.0;
  support_points[14](0) = vertices[6](0)/2.0+vertices[7](0)/2.0;
  support_points[14](1) = vertices[6](1)/2.0+vertices[7](1)/2.0;
  support_points[14](2) = vertices[6](2)/2.0+vertices[7](2)/2.0;
  support_points[15](0) = vertices[4](0)/2.0+vertices[7](0)/2.0;
  support_points[15](1) = vertices[4](1)/2.0+vertices[7](1)/2.0;
  support_points[15](2) = vertices[4](2)/2.0+vertices[7](2)/2.0;
  support_points[16](0) = vertices[0](0)/2.0+vertices[4](0)/2.0;
  support_points[16](1) = vertices[0](1)/2.0+vertices[4](1)/2.0;
  support_points[16](2) = vertices[0](2)/2.0+vertices[4](2)/2.0;
  support_points[17](0) = vertices[1](0)/2.0+vertices[5](0)/2.0;
  support_points[17](1) = vertices[1](1)/2.0+vertices[5](1)/2.0;
  support_points[17](2) = vertices[1](2)/2.0+vertices[5](2)/2.0;
  support_points[18](0) = vertices[2](0)/2.0+vertices[6](0)/2.0;
  support_points[18](1) = vertices[2](1)/2.0+vertices[6](1)/2.0;
  support_points[18](2) = vertices[2](2)/2.0+vertices[6](2)/2.0;
  support_points[19](0) = vertices[3](0)/2.0+vertices[7](0)/2.0;
  support_points[19](1) = vertices[3](1)/2.0+vertices[7](1)/2.0;
  support_points[19](2) = vertices[3](2)/2.0+vertices[7](2)/2.0;
  support_points[20](0) = vertices[0](0)/4.0+vertices[1](0)/4.0+vertices[2](0)/4.0+vertices[3](0)/4.0;
  support_points[20](1) = vertices[0](1)/4.0+vertices[1](1)/4.0+vertices[2](1)/4.0+vertices[3](1)/4.0;
  support_points[20](2) = vertices[0](2)/4.0+vertices[1](2)/4.0+vertices[2](2)/4.0+vertices[3](2)/4.0;
  support_points[21](0) = vertices[4](0)/4.0+vertices[5](0)/4.0+vertices[6](0)/4.0+vertices[7](0)/4.0;
  support_points[21](1) = vertices[4](1)/4.0+vertices[5](1)/4.0+vertices[6](1)/4.0+vertices[7](1)/4.0;
  support_points[21](2) = vertices[4](2)/4.0+vertices[5](2)/4.0+vertices[6](2)/4.0+vertices[7](2)/4.0;
  support_points[22](0) = vertices[0](0)/4.0+vertices[1](0)/4.0+vertices[4](0)/4.0+vertices[5](0)/4.0;
  support_points[22](1) = vertices[0](1)/4.0+vertices[1](1)/4.0+vertices[4](1)/4.0+vertices[5](1)/4.0;
  support_points[22](2) = vertices[0](2)/4.0+vertices[1](2)/4.0+vertices[4](2)/4.0+vertices[5](2)/4.0;
  support_points[23](0) = vertices[1](0)/4.0+vertices[2](0)/4.0+vertices[5](0)/4.0+vertices[6](0)/4.0;
  support_points[23](1) = vertices[1](1)/4.0+vertices[2](1)/4.0+vertices[5](1)/4.0+vertices[6](1)/4.0;
  support_points[23](2) = vertices[1](2)/4.0+vertices[2](2)/4.0+vertices[5](2)/4.0+vertices[6](2)/4.0;
  support_points[24](0) = vertices[2](0)/4.0+vertices[3](0)/4.0+vertices[6](0)/4.0+vertices[7](0)/4.0;
  support_points[24](1) = vertices[2](1)/4.0+vertices[3](1)/4.0+vertices[6](1)/4.0+vertices[7](1)/4.0;
  support_points[24](2) = vertices[2](2)/4.0+vertices[3](2)/4.0+vertices[6](2)/4.0+vertices[7](2)/4.0;
  support_points[25](0) = vertices[0](0)/4.0+vertices[3](0)/4.0+vertices[4](0)/4.0+vertices[7](0)/4.0;
  support_points[25](1) = vertices[0](1)/4.0+vertices[3](1)/4.0+vertices[4](1)/4.0+vertices[7](1)/4.0;
  support_points[25](2) = vertices[0](2)/4.0+vertices[3](2)/4.0+vertices[4](2)/4.0+vertices[7](2)/4.0;
  support_points[26](0) = vertices[0](0)/8.0+vertices[1](0)/8.0+vertices[2](0)/8.0+vertices[3](0)/8.0+vertices[4](0)/8.0+vertices[5](0)/8.0+vertices[6](0)/8.0+vertices[7](0)/8.0;
  support_points[26](1) = vertices[0](1)/8.0+vertices[1](1)/8.0+vertices[2](1)/8.0+vertices[3](1)/8.0+vertices[4](1)/8.0+vertices[5](1)/8.0+vertices[6](1)/8.0+vertices[7](1)/8.0;
  support_points[26](2) = vertices[0](2)/8.0+vertices[1](2)/8.0+vertices[2](2)/8.0+vertices[3](2)/8.0+vertices[4](2)/8.0+vertices[5](2)/8.0+vertices[6](2)/8.0+vertices[7](2)/8.0;
};



template <>
void FEQuadraticSub<3>::get_face_support_points (const typename DoFHandler<3>::face_iterator &face,
						 const Boundary<3>  &,
						 vector<Point<3> >  &support_points) const {
  Assert (support_points.size() == dofs_per_face,
	  ExcWrongFieldDimension (support_points.size(), dofs_per_face));

  for (unsigned int vertex=0; vertex<4; ++vertex)
    support_points[vertex] = face->vertex(vertex);
  for (unsigned int line=0; line<4; ++line)
    support_points[4+line] = (support_points[line] + support_points[(line+4)%4]) / 2;
  support_points[8] = (support_points[0] +
		       support_points[1] +
		       support_points[2] +
		       support_points[3]) / 4;
};

#endif





// explicit instantiations

template class FEQuadraticSub<deal_II_dimension>;

