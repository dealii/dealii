/* $Id$ */
/* Copyright W. Bangerth, University of Heidelberg, 1998 */


#include <fe/fe_lib.lagrange.h>
#include <grid/tria_iterator.h>
#include <grid/dof_accessor.h>
#include <grid/geometry_info.h>
#include <algorithm>




/*-----------------------------------------------------------------
 * For the 2D stuff, you may use the script below. However, apart
 * from the restriction matrices, I did not initially use it; it is
 * an adaption of the script for cubic and quartic elements. For
 * some of the data, however, a smaller script is given in the
 * FEQuadratic<2> constructor.
  n_functions      := 9:
  n_face_functions := 3:

  trial_function := (a1 + a2*xi + a3*xi*xi) +
                     (b1 + b2*xi + b3*xi*xi)*eta +
		     (c1 + c2*xi + c3*xi*xi)*eta*eta:
  face_trial_function := a + b*xi + c*xi*xi:
  # note: support_points[i] is a vector which is indexed from
  # one and not from zero!
  support_points := array(0..n_functions-1):
  support_points[0] := [0,0]:
  support_points[1] := [1,0]:
  support_points[2] := [1,1]:
  support_points[3] := [0,1]:
  support_points[4] := [1/2,0]:
  support_points[5] := [1,1/2]:
  support_points[6] := [1/2,1]:
  support_points[7] := [0,1/2]:
  support_points[8] := [1/2,1/2]:

  face_support_points := array(0..n_face_functions-1):
  face_support_points[0] := 0:
  face_support_points[1] := 1:
  face_support_points[2] := 1/2:

  constrained_face_support_points := array(0..2*(n_face_functions-2)+1-1):
  constrained_face_support_points[0] := 1/2:
  constrained_face_support_points[1] := 1/4:
  constrained_face_support_points[2] := 3/4:
  
  phi_polynom := array(0..n_functions-1):
  grad_phi_polynom := array(0..n_functions-1,0..1):
  local_mass_matrix := array(0..n_functions-1, 0..n_functions-1):
  prolongation := array(0..3,0..n_functions-1, 0..n_functions-1):
  interface_constraints := array(0..2*(n_face_functions-2)+1-1,
                                 0..n_face_functions-1):
  real_points := array(0..n_functions-1, 0..1);

  print ("Computing basis functions"):
  for i from 0 to n_functions-1 do
    print (i):
    values := array(1..n_functions):
    for j from 1 to n_functions do
      values[j] := 0:
    od:  
    values[i+1] := 1:

    equation_system := {}:
    for j from 0 to n_functions-1 do
      poly := subs(xi=support_points[j][1],
                   eta=support_points[j][2],
		   trial_function):
      if (i=j) then
        equation_system := equation_system union {poly = 1}:
      else	
        equation_system := equation_system union {poly = 0}:
      fi:	
    od:
    
    phi_polynom[i] := subs(solve(equation_system), trial_function):
    grad_phi_polynom[i,0] := diff(phi_polynom[i], xi):
    grad_phi_polynom[i,1] := diff(phi_polynom[i], eta):
  od:

  phi:= proc(i,x,y) subs(xi=x, eta=y, phi_polynom[i]): end:


  #points on children: let them be indexed one-based, as are
  #the support_points
  points[0] := array(0..n_functions-1, 1..2):
  points[1] := array(0..n_functions-1, 1..2):
  points[2] := array(0..n_functions-1, 1..2):
  points[3] := array(0..n_functions-1, 1..2):
  for i from 0 to n_functions-1 do
    points[0][i,1] := support_points[i][1]/2:
    points[0][i,2] := support_points[i][2]/2:
    
    points[1][i,1] := support_points[i][1]/2+1/2:
    points[1][i,2] := support_points[i][2]/2:

    points[2][i,1] := support_points[i][1]/2+1/2:
    points[2][i,2] := support_points[i][2]/2+1/2:

    points[3][i,1] := support_points[i][1]/2:
    points[3][i,2] := support_points[i][2]/2+1/2:
  od:  

  print ("Computing prolongation matrices"):
  for i from 0 to 3 do
    print ("child", i):
    for j from 0 to n_functions-1 do
      for k from 0 to n_functions-1 do
        prolongation[i,j,k] := phi(k, points[i][j,1], points[i][j,2]):
      od:
    od:
  od:

  print ("Computing restriction matrices"):
  # to get the restriction (interpolation) matrices, evaluate
  # the basis functions on the child cells at the global
  # interpolation points
  child_phi[0] := proc(i, x, y)
                    if ((x>1/2) or (y>1/2)) then
		      0:
		    else
		      phi(i,2*x,2*y):
		    fi:
		  end: 
  child_phi[1] := proc(i, x, y)
                    if ((x<1/2) or (y>1/2)) then
		      0:
		    else
		      phi(i,2*x-1,2*y):
		    fi:
		  end: 
  child_phi[2] := proc(i, x, y)
                    if ((x<1/2) or (y<1/2)) then
		      0:
		    else
		      phi(i,2*x-1,2*y-1):
		    fi:
		  end: 
  child_phi[3] := proc(i, x, y)
                    if ((x>1/2) or (y<1/2)) then
		      0:
		    else
		      phi(i,2*x,2*y-1):
		    fi:
		  end: 
  restriction := array(0..3,0..n_functions-1, 0..n_functions-1);  
  for child from 0 to 3 do
    for j from 0 to n_functions-1 do
      for k from 0 to n_functions-1 do
        restriction[child,j,k] := child_phi[child](k,
	                                           support_points[j][1],
						   support_points[j][2]):
      od:
    od:
  od:

  
  print ("Computing local mass matrix"):
  # tphi are the basis functions of the linear element. These functions
  # are used for the computation of the subparametric transformation from
  # unit cell to real cell.
  # x and y are arrays holding the x- and y-values of the four vertices
  # of this cell in real space. 
  x := array(0..3);
  y := array(0..3);
  tphi[0] := (1-xi)*(1-eta):
  tphi[1] := xi*(1-eta):
  tphi[2] := xi*eta:
  tphi[3] := (1-xi)*eta:
  x_real := sum(x[s]*tphi[s], s=0..3):
  y_real := sum(y[s]*tphi[s], s=0..3):
  detJ := diff(x_real,xi)*diff(y_real,eta) - diff(x_real,eta)*diff(y_real,xi):
  for i from 0 to n_functions-1 do
    print ("line", i):
    for j from 0 to n_functions-1 do
      local_mass_matrix[i,j] := int(int(phi_polynom[i] * phi_polynom[j] * detJ,
                                        xi=0..1), eta=0..1):
    od:
  od:

  print ("computing support points in real space"):
  for i from 0 to n_functions-1 do
    real_points[i,0] := subs(xi=support_points[i][1],
                             eta=support_points[i][2], x_real);
    real_points[i,1] := subs(xi=support_points[i][1],
                             eta=support_points[i][2], y_real);
  od:

  print ("computing interface constraint matrices"):
  # compute the interface constraint matrices. these are the values of the
  # basis functions on the coarse cell's face at the points of the child
  # cell's basis functions on the child faces
  face_phi_polynom := array(0..n_face_functions-1):
  for i from 0 to n_face_functions-1 do
    # note that the interp function wants vectors indexed from
    #   one and not from zero. 
    values := array(1..n_face_functions):
    for j from 1 to n_face_functions do
      values[j] := 0:
    od:  
    values[i+1] := 1:

    shifted_face_support_points := array (1..n_face_functions):
    for j from 1 to n_face_functions do
      shifted_face_support_points[j] := face_support_points[j-1]:
    od:
    
    face_phi_polynom[i] := interp (shifted_face_support_points, values, xi):
  od:

  for i from 0 to 2*(n_face_functions-2)+1-1 do
    for j from 0 to n_face_functions-1 do
      interface_constraints[i,j] := subs(xi=constrained_face_support_points[i],
                                     face_phi_polynom[j]); 
    od:
  od:


  print ("writing data to files"):
  readlib(C):
  C(phi_polynom, filename=shape_value_2d):
  C(grad_phi_polynom, filename=shape_grad_2d):
  C(prolongation, filename=prolongation_2d):
  C(restriction, filename=restriction_2d):
  C(local_mass_matrix, optimized, filename=massmatrix_2d):
  C(interface_constraints, filename=constraints_2d):
  C(real_points, optimized, filename=real_points_2d);

  ---------------------------------------------------------------

  Use the following perl scripts to convert the output into a
  suitable format.
  
  perl -pi -e 's/phi_polynom\[(\d+)\] =/case $1: return/g;' shape_value_2d
  perl -pi -e 's/([^;])\n/$1/g;' shape_grad_2d
  perl -pi -e 's/grad_phi_polynom\[(\d+)\]\[0\] = (.*);/case $1: return Point<2>($2,/g;' shape_grad_2d
  perl -pi -e 's/grad_phi_polynom\[(\d+)\]\[1\] = (.*);/$2);/g;' shape_grad_2d
  perl -pi -e 's/\[(\d+)\]\[(\d+)\]/($1,$2)/g;' massmatrix_2d
  perl -pi -e 's/(t\d+) =/const double $1 =/g;' massmatrix_2d
  perl -pi -e 's/\[(\d+)\]\[(\d+)\]\[(\d+)\]/[$1]($2,$3)/g;' prolongation_2d
  perl -pi -e 's/.*= 0.0;\n//g;' prolongation_2d
  perl -pi -e 's/\[(\d+)\]\[(\d+)\]\[(\d+)\]/[$1]($2,$3)/g;' restriction_2d
  perl -pi -e 's/.*= 0.0;\n//g;' restriction_2d
  perl -pi -e 's/\[(\d+)\]\[(\d+)\]/($1,$2)/g;' constraints_2d
*/





#if deal_II_dimension == 1

template <>
FEQuadraticSub<1>::FEQuadraticSub () :
		FELinearMapping<1> (1, 1) {
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
Point<1>
FEQuadraticSub<1>::shape_grad(const unsigned int i,
			      const Point<1>    &p) const
{
  Assert((i<total_dofs), ExcInvalidIndex(i));
  const double xi = p(0);
  switch (i)
    {
      case 0: return Point<1>(-3+4*xi);
      case 1: return Point<1>(4*xi-1);
      case 2: return Point<1>(4-8*xi);
    }
  return Point<1>();
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
Point<2>
FEQuadraticSub<2>::shape_grad (const unsigned int i,
			       const Point<2>    &p) const
{
  Assert (i<total_dofs, ExcInvalidIndex(i));

  const double xi = p(0),
	       eta= p(1);
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





// explicit instantiations

template class FEQuadraticSub<deal_II_dimension>;

