/* $Id$ */

#include <fe/fe_lib.h>
#include <grid/tria_iterator.h>
#include <grid/dof_accessor.h>
#include <grid/geometry_info.h>
#include <algorithm>


/*--------------------------------- For 1d ---------------------------------
  -- Use the following maple script to generate the basis functions,
  -- gradients and prolongation matrices as well as the mass matrix.
  -- Make sure that the files do not exists beforehand, since output
  -- is appended instead of overwriting previous contents.
  --
  -- You should only have to change the very first lines for polynomials
  -- of higher order.
  --------------------------------------------------------------------------
  n_functions := 4;
  
  ansatz_points := array(0..n_functions-1);
  ansatz_points[0] := 0;
  ansatz_points[1] := 1;
  ansatz_points[2] := 1/3;
  ansatz_points[3] := 2/3;

  phi_polynom := array(0..n_functions-1);
  grad_phi_polynom := array(0..n_functions-1);
  local_mass_matrix := array(0..n_functions-1, 0..n_functions-1);

  for i from 0 to n_functions-1 do
    # note that the interp function wants vector indexed from
    #   one and not from zero. 
    values := array(1..n_functions);
    for j from 1 to n_functions do
      values[j] := 0;
    od;  
    values[i+1] := 1;

    shifted_ansatz_points := array (1..n_functions);
    for j from 1 to n_functions do
      shifted_ansatz_points[j] := ansatz_points[j-1];
    od;
    
    phi_polynom[i] := interp (shifted_ansatz_points, values, xi);
    grad_phi_polynom[i] := diff(phi_polynom[i], xi);
  od;

  phi:= proc(i,x) subs(xi=x, phi_polynom[i]); end;


  points[0] := array(0..n_functions-1);
  points[1] := array(0..n_functions-1);
  for i from 0 to n_functions-1 do
    points[0][i] := ansatz_points[i]/2;  
    points[1][i] := ansatz_points[i]/2+1/2;
  od;  

  prolongation := array(0..1,0..n_functions-1, 0..n_functions-1);

  for i from 0 to 1 do
    for j from 0 to n_functions-1 do
      for k from 0 to n_functions-1 do
        prolongation[i,j,k] := phi(k, points[i][j]);
      od;
    od;
  od;

  for i from 0 to n_functions-1 do
    for j from 0 to n_functions-1 do
      local_mass_matrix[i,j] := int(phi_polynom[i] * phi_polynom[j] * h,
                                    xi=0..1);
    od;
  od;
  
  readlib(C);
  C(phi_polynom, filename=shape_value_1d);
  C(grad_phi_polynom, filename=shape_grad_1d);
  C(prolongation, filename=prolongation_1d);
  C(local_mass_matrix, optimized, filename=massmatrix_1d);

  -----------------------------------------------------------------------
  Use the following perl scripts to convert the output into a
  suitable format:
  
  perl -pi -e 's/phi_polynom\[(\d)\] =/case $1: return/g;' shape_value_1d
  perl -pi -e 's/grad_phi_polynom\[(\d)\] = (.*);/case $1: return Point<1>($2);/g;' shape_grad_1d
  perl -pi -e 's/\[(\d)\]\[(\d)\]/($1,$2)/g;' massmatrix_1d
  perl -pi -e 's/\[(\d)\]\[(\d)\]\[(\d)\]/[$1]($2,$3)/g;' prolongation_1d
  perl -pi~ -e 's/(t\d)/const double $1/g;' massmatrix_1d
*/




/*--------------------------------- For 2d ---------------------------------
  -- Use the following maple script to generate the basis functions,
  -- gradients and prolongation matrices as well as the mass matrix.
  -- Make sure that the files do not exists beforehand, since output
  -- is appended instead of overwriting previous contents.
  --
  -- You should only have to change the very first lines for polynomials
  -- of higher order.
  --------------------------------------------------------------------------
  n_functions := 16;

  ansatz_function := (a1 + a2*xi + a3*xi*xi + a4*xi*xi*xi) +
                     (b1 + b2*xi + b3*xi*xi + b4*xi*xi*xi)*eta +
		     (c1 + c2*xi + c3*xi*xi + c4*xi*xi*xi)*eta*eta +
		     (d1 + d2*xi + d3*xi*xi + d4*xi*xi*xi)*eta*eta*eta;
  ansatz_points := array(0..n_functions-1);
  # note: ansatz_points[i] is a vector which is indexed from
  # one and not from zero!
  ansatz_points[0] := [0,0];
  ansatz_points[1] := [1,0];
  ansatz_points[2] := [1,1];
  ansatz_points[3] := [0,1];
  ansatz_points[4] := [1/3,0];
  ansatz_points[5] := [2/3,0];
  ansatz_points[6] := [1,1/3];
  ansatz_points[7] := [1,2/3];
  ansatz_points[8] := [1/3,1];
  ansatz_points[9] := [2/3,1];
  ansatz_points[10]:= [0,1/3];
  ansatz_points[11]:= [0,2/3];
  ansatz_points[12]:= [1/3,1/3];
  ansatz_points[13]:= [2/3,1/3];
  ansatz_points[14]:= [2/3,2/3];
  ansatz_points[15]:= [1/3,2/3];

  
  phi_polynom := array(0..n_functions-1);
  grad_phi_polynom := array(0..n_functions-1,0..1);
  local_mass_matrix := array(0..n_functions-1, 0..n_functions-1);
  prolongation := array(0..3,0..n_functions-1, 0..n_functions-1);


  for i from 0 to n_functions-1 do
    values := array(1..n_functions);
    for j from 1 to n_functions do
      values[j] := 0;
    od;  
    values[i+1] := 1;

    equation_system := {};
    for j from 0 to n_functions-1 do
      poly := subs(xi=ansatz_points[j][1],
                   eta=ansatz_points[j][2],
		   ansatz_function);
      if (i=j) then
        equation_system := equation_system union {poly = 1};
      else	
        equation_system := equation_system union {poly = 0};
      fi;	
    od;
    
    phi_polynom[i] := subs(solve(equation_system), ansatz_function);
    grad_phi_polynom[i,0] := diff(phi_polynom[i], xi);
    grad_phi_polynom[i,1] := diff(phi_polynom[i], eta);
  od;

  phi:= proc(i,x,y) subs(xi=x, eta=y, phi_polynom[i]); end;

  #points on children; let them be indexed one-based, as are
  #the ansatz_points
  points[0] := array(0..n_functions-1, 1..2);
  points[1] := array(0..n_functions-1, 1..2);
  points[2] := array(0..n_functions-1, 1..2);
  points[3] := array(0..n_functions-1, 1..2);
  for i from 0 to n_functions-1 do
    points[0][i,1] := ansatz_points[i][1]/2;
    points[0][i,2] := ansatz_points[i][2]/2;
    
    points[1][i,1] := ansatz_points[i][1]/2+1/2;
    points[1][i,2] := ansatz_points[i][2]/2;

    points[2][i,1] := ansatz_points[i][1]/2+1/2;
    points[2][i,2] := ansatz_points[i][2]/2+1/2;

    points[3][i,1] := ansatz_points[i][1]/2;
    points[3][i,2] := ansatz_points[i][2]/2+1/2;
  od;  

  for i from 0 to 3 do
    for j from 0 to n_functions-1 do
      for k from 0 to n_functions-1 do
        prolongation[i,j,k] := phi(k, points[i][j,1], points[i][j,2]);
      od;
    od;
  od;

  # tphi are the basis functions of the linear element. These functions
  # are used for the computation of the subparametric transformation from
  # unit cell to real cell.
  tphi[0] := (1-xi)*(1-eta);
  tphi[1] := xi*(1-eta);
  tphi[2] := xi*eta;
  tphi[3] := (1-xi)*eta;
  x_real := sum(x[s]*tphi[s], s=0..3);
  y_real := sum(y[s]*tphi[s], s=0..3);
  detJ := diff(x_real,xi)*diff(y_real,eta) - diff(x_real,eta)*diff(y_real,xi);
  for i from 0 to n_functions-1 do
    for j from 0 to n_functions-1 do
      local_mass_matrix[i,j] := int(int(phi_polynom[i] * phi_polynom[j] * detJ,
                                        xi=0..1), eta=0..1);
    od;
  od;
  
  readlib(C);
  C(phi_polynom, filename=shape_value_2d);
  C(grad_phi_polynom, filename=shape_grad_2d);
  C(prolongation, filename=prolongation_2d);
  C(local_mass_matrix, optimized, filename=massmatrix_2d);

  -----------------------------------------------------------------------
  Use the following perl scripts to convert the output into a
  suitable format.
  
  perl -pi -e 's/phi_polynom\[(\d+)\] =/case $1: return/g;' shape_value_2d
  perl -pi -e 's/([^;])\n/$1/g;' shape_grad_2d
  perl -pi -e 's/grad_phi_polynom\[(\d+)\]\[0\] = (.*);/case $1: return Point<2>($2,/g;' shape_grad_2d
  perl -pi -e 's/grad_phi_polynom\[(\d+)\]\[1\] = (.*);/$2);/g;' shape_grad_2d
  perl -pi -e 's/\[(\d+)\]\[(\d+)\]/($1,$2)/g;' massmatrix_2d
  perl -pi -e 's/\[(\d+)\]\[(\d+)\]\[(\d+)\]/[$1]($2,$3)/g;' prolongation_2d
  perl -pi~ -e 's/(t\d+) =/const double $1 =/g;' massmatrix_2d
*/






#if deal_II_dimension == 1

template <>
FECubicSub<1>::FECubicSub () :
		FiniteElement<1> (1, 2) {
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
};



template <>
void FECubicSub<1>::fill_fe_values (const DoFHandler<1>::cell_iterator &cell,
				     const vector<Point<1> >            &unit_points,
				     vector<dFMatrix>  &jacobians,
				     const bool         compute_jacobians,
				     vector<Point<1> > &ansatz_points,
				     const bool         compute_ansatz_points,
				     vector<Point<1> > &q_points,
				     const bool         compute_q_points,
				     const Boundary<1> &boundary) const {
				   // simply pass down
  FiniteElement<1>::fill_fe_values (cell, unit_points,
				    jacobians, compute_jacobians,
				    ansatz_points, compute_ansatz_points,
				    q_points, compute_q_points, boundary);
};



template <>
double
FECubicSub<1>::shape_value(const unsigned int i,
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
inline
double
FECubicSub<1>::linear_shape_value(const unsigned int i,
				      const Point<1>     &p) const
{
  Assert((i<2), ExcInvalidIndex(i));
  const double xi = p(0);
  switch (i)
    {
      case 0: return 1.-xi;
      case 1: return xi;
    }
  return 0.;
};



template <>
Point<1>
FECubicSub<1>::shape_grad(const unsigned int i,
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
inline
Point<1>
FECubicSub<1>::linear_shape_grad(const unsigned int i,
				     const Point<1>&) const
{
  Assert((i<2), ExcInvalidIndex(i));
  switch (i)
    {
    case 0: return Point<1>(-1.);
    case 1: return Point<1>(1.);
    }
  return Point<1>();
};



template <>
void FECubicSub<1>::get_ansatz_points (const typename DoFHandler<1>::cell_iterator &cell,
					   const Boundary<1>  &boundary,
					   vector<Point<1> >  &ansatz_points) const {
  FiniteElement<1>::get_ansatz_points (cell, boundary, ansatz_points);
};



template <>
void FECubicSub<1>::get_face_ansatz_points (const typename DoFHandler<1>::face_iterator &,
					     const Boundary<1>  &,
					     vector<Point<1> >  &) const {
  Assert (false, ExcInternalError());
};



template <>
void FECubicSub<1>::get_face_jacobians (const DoFHandler<1>::face_iterator &,
					    const Boundary<1>         &,
					    const vector<Point<0> > &,
					    vector<double>      &) const {
  Assert (false, ExcInternalError());
};



template <>
void FECubicSub<1>::get_subface_jacobians (const DoFHandler<1>::face_iterator &,
					       const unsigned int           ,
					       const vector<Point<0> > &,
					       vector<double>      &) const {
  Assert (false, ExcInternalError());
};



template <>
void FECubicSub<1>::get_normal_vectors (const DoFHandler<1>::cell_iterator &,
					    const unsigned int,
					    const Boundary<1> &,
					    const vector<Point<0> > &,
					    vector<Point<1> > &) const {
  Assert (false, ExcInternalError());
};



template <>
void FECubicSub<1>::get_normal_vectors (const DoFHandler<1>::cell_iterator &,
					    const unsigned int,
					    const unsigned int,
					    const vector<Point<0> > &,
					    vector<Point<1> > &) const {
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
  local_mass_matrix(0,0) =  const double t1;
  local_mass_matrix(0,1) =  const double t2;
  local_mass_matrix(0,2) =  const double t3;
  local_mass_matrix(0,3) = - const double t4;
  local_mass_matrix(1,0) =  const double t2;
  local_mass_matrix(1,1) =  const double t1;
  local_mass_matrix(1,2) = - const double t4;
  local_mass_matrix(1,3) =  const double t3;
  local_mass_matrix(2,0) =  const double t3;
  local_mass_matrix(2,1) = - const double t4;
  local_mass_matrix(2,2) =  const double t5;
  local_mass_matrix(2,3) = - const double t6;
  local_mass_matrix(3,0) = - const double t4;
  local_mass_matrix(3,1) =  const double t3;
  local_mass_matrix(3,2) = - const double t6;
  local_mass_matrix(3,3) =  const double t5;
};

#endif




#if deal_II_dimension == 2

template <>
FECubicSub<2>::FECubicSub () :
		FiniteElement<2> (1, 1, 1)
{
  interface_constraints(0,2) = 1.0;
  interface_constraints(1,0) = 3./8.;
  interface_constraints(1,1) = -1./8.;
  interface_constraints(1,2) = 3./4.;
  interface_constraints(2,0) = -1./8.;
  interface_constraints(2,1) = 3./8.;
  interface_constraints(2,2) = 3./4.;

/*
  Get the prolongation matrices by the following little maple script:

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
inline
double
FECubicSub<2>::linear_shape_value (const unsigned int i,
				       const Point<2>& p) const
{
  Assert((i<4), ExcInvalidIndex(i));
  switch (i)
    {
    case 0: return (1.-p(0)) * (1.-p(1));
    case 1: return p(0) * (1.-p(1));
    case 2: return p(0) * p(1);
    case 3: return (1.-p(0)) * p(1);
    }
  return 0.;
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
inline
Point<2>
FECubicSub<2>::linear_shape_grad (const unsigned int i,
				      const Point<2>& p) const
{
  Assert((i<4), ExcInvalidIndex(i));
  switch (i)
    {
    case 0: return Point<2> (p(1)-1., p(0)-1.);
    case 1: return Point<2> (1.-p(1), -p(0));
    case 2: return Point<2> (p(1), p(0));
    case 3: return Point<2> (-p(1), 1.-p(0));
    }
  return Point<2> ();
};



template <>
void FECubicSub<2>::get_local_mass_matrix (const DoFHandler<2>::cell_iterator &cell,
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
void FECubicSub<2>::get_ansatz_points (const typename DoFHandler<2>::cell_iterator &cell,
					   const Boundary<2>&,
					   vector<Point<2> >  &ansatz_points) const {
  Assert (ansatz_points.size() == total_dofs,
	  ExcWrongFieldDimension (ansatz_points.size(), total_dofs));
  
  for (unsigned int vertex=0; vertex<4; ++vertex)
    ansatz_points[vertex] = cell->vertex(vertex);

				   // for the bilinear mapping, the centers
				   // of the face on the unit cell are mapped
				   // to the mean coordinates of the vertices
  for (unsigned int line=0; line<4; ++line)
    ansatz_points[4+line] = (cell->line(line)->vertex(0) +
			     cell->line(line)->vertex(1)) / 2;
				   // same for the center of the square:
				   // since all four linear basis functions
				   // take on the value 1/4 at the center,
				   // the center is mapped to the mean
				   // coordinates of the four vertices
  ansatz_points[8] = (ansatz_points[0] +
		      ansatz_points[1] +
		      ansatz_points[2] +
		      ansatz_points[3]) / 4;
};



template <>
void FECubicSub<2>::get_face_ansatz_points (const typename DoFHandler<2>::face_iterator &face,
						const Boundary<2>  &,
						vector<Point<2> >  &ansatz_points) const {
  Assert (ansatz_points.size() == dofs_per_face,
	  ExcWrongFieldDimension (ansatz_points.size(), dofs_per_face));

  for (unsigned int vertex=0; vertex<2; ++vertex)
    ansatz_points[vertex] = face->vertex(vertex);
  ansatz_points[2] = (ansatz_points[0] + ansatz_points[1]) / 2;
};



template <>
void FECubicSub<2>::get_face_jacobians (const DoFHandler<2>::face_iterator &face,
					    const Boundary<2>         &,
					    const vector<Point<1> > &unit_points,
					    vector<double> &face_jacobians) const {
				   // more or less copied from the linear
				   // finite element
  Assert (unit_points.size() == face_jacobians.size(),
	  ExcWrongFieldDimension (unit_points.size(), face_jacobians.size()));

				   // a linear mapping for a single line
				   // produces particularly simple
				   // expressions for the jacobi
				   // determinant :-)
  const double h = sqrt((face->vertex(1) - face->vertex(0)).square());
  fill_n (face_jacobians.begin(),
	  unit_points.size(),
	  h);  
};



template <>
void FECubicSub<2>::get_subface_jacobians (const DoFHandler<2>::face_iterator &face,
					      const unsigned int           ,
					      const vector<Point<1> > &unit_points,
					      vector<double> &face_jacobians) const {
				   // more or less copied from the linear
				   // finite element
  Assert (unit_points.size() == face_jacobians.size(),
	  ExcWrongFieldDimension (unit_points.size(), face_jacobians.size()));
  Assert (face->at_boundary() == false,
	  ExcBoundaryFaceUsed ());

				   // a linear mapping for a single line
				   // produces particularly simple
				   // expressions for the jacobi
				   // determinant :-)
  const double h = sqrt((face->vertex(1) - face->vertex(0)).square());
  fill_n (face_jacobians.begin(),
	  unit_points.size(),
	  h/2);
};



template <>
void FECubicSub<2>::get_normal_vectors (const DoFHandler<2>::cell_iterator &cell,
					   const unsigned int       face_no,
					   const Boundary<2>       &,
					   const vector<Point<1> > &unit_points,
					   vector<Point<2> > &normal_vectors) const {
				   // more or less copied from the linear
				   // finite element
  Assert (unit_points.size() == normal_vectors.size(),
	  ExcWrongFieldDimension (unit_points.size(), normal_vectors.size()));

  const DoFHandler<2>::face_iterator face = cell->face(face_no);
				   // compute direction of line
  const Point<2> line_direction = (face->vertex(1) - face->vertex(0));
				   // rotate to the right by 90 degrees
  const Point<2> normal_direction(line_direction(1),
				  -line_direction(0));

  if (face_no <= 1)
				     // for sides 0 and 1: return the correctly
				     // scaled vector
    fill (normal_vectors.begin(), normal_vectors.end(),
	  normal_direction / sqrt(normal_direction.square()));
  else
				     // for sides 2 and 3: scale and invert
				     // vector
    fill (normal_vectors.begin(), normal_vectors.end(),
	  normal_direction / (-sqrt(normal_direction.square())));
};



template <>
void FECubicSub<2>::get_normal_vectors (const DoFHandler<2>::cell_iterator &cell,
					   const unsigned int       face_no,
					   const unsigned int,
					   const vector<Point<1> > &unit_points,
					   vector<Point<2> > &normal_vectors) const {
				   // more or less copied from the linear
				   // finite element
				   // note, that in 2D the normal vectors to the
				   // subface have the same direction as that
				   // for the face
  Assert (unit_points.size() == normal_vectors.size(),
	  ExcWrongFieldDimension (unit_points.size(), normal_vectors.size()));
  Assert (cell->face(face_no)->at_boundary() == false,
	  ExcBoundaryFaceUsed ());

  const DoFHandler<2>::face_iterator face = cell->face(face_no);
				   // compute direction of line
  const Point<2> line_direction = (face->vertex(1) - face->vertex(0));
				   // rotate to the right by 90 degrees
  const Point<2> normal_direction(line_direction(1),
				  -line_direction(0));

  if (face_no <= 1)
				     // for sides 0 and 1: return the correctly
				     // scaled vector
    fill (normal_vectors.begin(), normal_vectors.end(),
	  normal_direction / sqrt(normal_direction.square()));
  else
				     // for sides 2 and 3: scale and invert
				     // vector
    fill (normal_vectors.begin(), normal_vectors.end(),
	  normal_direction / (-sqrt(normal_direction.square())));
};

#endif





template <int dim>
void FECubicSub<dim>::fill_fe_values (const DoFHandler<dim>::cell_iterator &cell,
					  const vector<Point<dim> >            &unit_points,
					  vector<dFMatrix>    &jacobians,
					  const bool           compute_jacobians,
					  vector<Point<dim> > &ansatz_points,
					  const bool           compute_ansatz_points,
					  vector<Point<dim> > &q_points,
					  const bool           compute_q_points,
					  const Boundary<dim> &boundary) const {
  Assert (jacobians.size() == unit_points.size(),
	  ExcWrongFieldDimension(jacobians.size(), unit_points.size()));
  Assert (q_points.size() == unit_points.size(),
	  ExcWrongFieldDimension(q_points.size(), unit_points.size()));
  Assert (ansatz_points.size() == total_dofs,
	  ExcWrongFieldDimension(ansatz_points.size(), total_dofs));

  
  unsigned int n_points=unit_points.size();

  Point<dim> vertices[GeometryInfo<dim>::vertices_per_cell];
  for (unsigned int l=0; l<GeometryInfo<dim>::vertices_per_cell; ++l)
    vertices[l] = cell->vertex(l);
  

  if (compute_q_points) 
    {
				       // initialize points to zero
      for (unsigned int i=0; i<n_points; ++i)
	q_points[i] = Point<dim> ();
      
				       // note: let x_l be the vector of the
				       // lth quadrature point in real space and
				       // xi_l that on the unit cell, let further
				       // p_j be the vector of the jth vertex
				       // of the cell in real space and
				       // N_j(xi_l) be the value of the associated
				       // basis function at xi_l, then
				       // x_l(xi_l) = sum_j p_j N_j(xi_l)
				       //
				       // Here, N_j is the *linear* basis function,
				       // not that of the finite element, since we
				       // use a subparametric mapping
      for (unsigned int j=0; j<GeometryInfo<dim>::vertices_per_cell; ++j) 
	for (unsigned int l=0; l<n_points; ++l) 
	  q_points[l] += vertices[j] * linear_shape_value(j, unit_points[l]);
    };
  

/* jacobi matrices: compute d(x)/d(xi) and invert this
   Let M(l) be the inverse of J at the quadrature point l, then
     M_{ij}(l) = sum_s p_i(s) d(N_s(l))/d(xi_j)
   where p_i(s) is the i-th coordinate of the s-th vertex vector,
   N_s(l) is the value of the s-th vertex shape function at the
   quadrature point l.

   We could therefore write:
   l=0..n_points-1
     i=0..dim-1
       j=0..dim-1
         M_{ij}(l) = 0
	 s=0..n_vertices
	   M_{ij}(l) += p_i(s) d(N_s(l))/d(xi_j)

  However, we rewrite the loops to only compute the gradient once for
  each integration point and basis function.
*/
  if (compute_jacobians) 
    {
      dFMatrix M(dim,dim);
      for (unsigned int l=0; l<n_points; ++l) 
	{
	  M.clear ();
	  for (unsigned int s=0; s<GeometryInfo<dim>::vertices_per_cell; ++s)
	    {
					       // we want the linear transform,
					       // so use that function
	      const Point<dim> gradient = linear_shape_grad (s, unit_points[l]);
	      for (unsigned int i=0; i<dim; ++i)
		for (unsigned int j=0; j<dim; ++j)
		  M(i,j) += vertices[s](i) * gradient(j);
	    };
	  jacobians[l].invert(M);
	};
    };

				   // compute ansatz points, which are
				   // the corners for linear elements
  if (compute_ansatz_points)
    get_ansatz_points (cell, boundary, ansatz_points);
};




// explicit instantiations

template class FECubicSub<deal_II_dimension>;

