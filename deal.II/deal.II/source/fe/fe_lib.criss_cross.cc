/* $Id$ */
/* Copyright W. Bangerth, University of Heidelberg, 1998 */


#include <fe/fe_lib.criss_cross.h>
#include <grid/tria_iterator.h>
#include <grid/dof_accessor.h>
#include <grid/geometry_info.h>




/*-----------------------------------2d------------------------------------
  Maple script to automate some of the error-prone computations on
  this strange sort of element.

  n_functions      := 5:

  # note: support_points[i] is a vector which is indexed from
  # one and not from zero!
  support_points[0] := [0,0]:
  support_points[1] := [1,0]:
  support_points[2] := [1,1]:
  support_points[3] := [0,1]:
  support_points[4] := [1/2,1/2]:
  
  phi[0] := proc(x,y) if(y<1-x) then 1-x-y; else 0; fi; end:
  phi[1] := proc(x,y) if(y<x)   then x-y;   else 0; fi; end:
  phi[2] := proc(x,y) if(y>1-x) then x+y-1; else 0; fi; end:
  phi[3] := proc(x,y) if(y>x)   then y-x;   else 0; fi; end:
  phi[4] := proc(x,y) 1 - phi[0](x,y) - phi[1](x,y)
                        - phi[2](x,y) - phi[3](x,y) ; end:

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

  prolongation := array(0..3,0..n_functions-1, 0..n_functions-1):
  print ("Computing prolongation matrices"):
  for i from 0 to 3 do
    print ("child", i):
    for j from 0 to n_functions-1 do
      for k from 0 to n_functions-1 do
        prolongation[i,j,k] := phi[k](points[i][j,1], points[i][j,2]):
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
		      phi[i](2*x,2*y):
		    fi:
		  end: 
  child_phi[1] := proc(i, x, y)
                    if ((x<1/2) or (y>1/2)) then
		      0:
		    else
		      phi[i](2*x-1,2*y):
		    fi:
		  end: 
  child_phi[2] := proc(i, x, y)
                    if ((x<1/2) or (y<1/2)) then
		      0:
		    else
		      phi[i](2*x-1,2*y-1):
		    fi:
		  end: 
  child_phi[3] := proc(i, x, y)
                    if ((x>1/2) or (y<1/2)) then
		      0:
		    else
		      phi[i](2*x,2*y-1):
		    fi:
		  end: 
  restriction := array(0..3,0..n_functions-1, 0..n_functions-1):
  for child from 0 to 3 do
    for j from 0 to n_functions-1 do
      for k from 0 to n_functions-1 do
        restriction[child,j,k] := child_phi[child](k,
	                                           support_points[j][1],
						   support_points[j][2]):
      od:
    od:
  od:

  
  # these are the basis functions differentiated with respect to
  # xi and eta. we need them for the computation of the jacobi
  # matrix, since we can't just differentiate a function.
  phi_xi[0] := proc(x,y) if(y<1-x) then -1;  else 0; fi; end:
  phi_xi[1] := proc(x,y) if(y<x)   then 1;   else 0; fi; end:
  phi_xi[2] := proc(x,y) if(y>1-x) then 1;   else 0; fi; end:
  phi_xi[3] := proc(x,y) if(y>x)   then -1;  else 0; fi; end:
  phi_xi[4] := proc(x,y) 1 - phi_xi[0](x,y) - phi_xi[1](x,y)
                           - phi_xi[2](x,y) - phi_xi[3](x,y) ; end:

  phi_eta[0] := proc(x,y) if(y<1-x) then -1;  else 0; fi; end:
  phi_eta[1] := proc(x,y) if(y<x)   then -1;  else 0; fi; end:
  phi_eta[2] := proc(x,y) if(y>1-x) then 1;   else 0; fi; end:
  phi_eta[3] := proc(x,y) if(y>x)   then 1;   else 0; fi; end:
  phi_eta[4] := proc(x,y) 1 - phi_eta[0](x,y) - phi_eta[1](x,y)
                            - phi_eta[2](x,y) - phi_eta[3](x,y) ; end:

  # define an array of the support points in real space; the first
  # four are the vertices, the last one is the crossing point of
  # the two diagonals
  print ("Computing cross point"):
  x := array(0..4):
  y := array(0..4):

  eq_sys := {(1-t)*x[0] + t*x[2] = (1-s)*x[1] + s*x[3],
             (1-t)*y[0] + t*y[2] = (1-s)*y[1] + s*y[3]}:
  solution := solve (eq_sys, {s,t}):

  # set last point in dependence of the first four
  x[4] := subs (solution, (1-t)*x[0] + t*x[2]):
  y[4] := subs (solution, (1-t)*y[0] + t*y[2]):

  # this is the mapping from the unit cell to the real cell, only for
  # completeness; we can't use it here, since phi[i] can't be
  # differentiated
  x_real := simplify(sum(x[s]*phi[s], s=0..4)):
  y_real := simplify(sum(y[s]*phi[s], s=0..4)):

  # correct form of the jacobi determinant:
  #   detJ :=   diff(x_real,xi)*diff(y_real,eta)
  #           - diff(x_real,eta)*diff(y_real,xi):
  # better now:
  detJ1 := proc(xi,eta) sum(x[s]*phi_xi[s](xi,eta), s=0..4); end:
  detJ2 := proc(xi,eta) sum(y[s]*phi_eta[s](xi,eta), s=0..4); end:
  detJ3 := proc(xi,eta) sum(x[s]*phi_eta[s](xi,eta), s=0..4); end:
  detJ4 := proc(xi,eta) sum(y[s]*phi_xi[s](xi,eta), s=0..4); end:
  detJ := proc(xi,eta)
             detJ1(xi,eta) * detJ2(xi,eta) -
	     detJ3(xi,eta) * detJ4(xi,eta);
          end:


  # Now for the mass matrix: we divide the entire cell into four
  # sectors:
  #
  # *-------*
  # |\     /|
  # | \ 3 / |
  # |  \ /  |
  # |4  *  2|
  # |  / \  |
  # | / 1 \ |
  # |/     \|
  # *-------*
  #
  # In each of these sectors, the Jacobi determinant is constant
  # so that we can assemble the local mass matrix by summation
  # over these four sectors. Since the basis functions are as of
  # now only expressed as if-then-else statements, we have to
  # express them for each sector separately and name them
  # phi_s[i]. detJ_s denotes the Jacobi determinant on this sector.

  print ("Computing mass matrix"):

  mass_matrix := array (0..n_functions-1, 0..n_functions-1):
  for i from 0 to n_functions-1 do
    for j from 0 to n_functions-1 do
      mass_matrix[i,j] := 0:
    od:
  od:
  
  # sector 1
  phi_s[0] := 1-x-y:
  phi_s[1] := x-y:
  phi_s[2] := 0:
  phi_s[3] := 0:
  phi_s[4] := 1 - phi_s[0] - phi_s[1] - phi_s[2] - phi_s[3]:

  detJ_s := simplify(detJ(1/2, 1/4)):
  
  for i from 0 to n_functions-1 do
    for j from 0 to n_functions-1 do
      # split integral over sector into the two parts
      mass_matrix[i,j] := mass_matrix[i,j] +
                          int(int(phi_s[i] * phi_s[j],
                                  y=0..x),
			      x=0..1/2) * detJ_s:
      mass_matrix[i,j] := mass_matrix[i,j] +
                          int(int(phi_s[i] * phi_s[j],
                                  y=0..1-x),
			      x=1/2..1) * detJ_s:
    od:
  od:  
  
  # sector 2
  phi_s[0] := 0:
  phi_s[1] := x-y:
  phi_s[2] := x+y-1:
  phi_s[3] := 0:
  phi_s[4] := 1 - phi_s[0] - phi_s[1] - phi_s[2] - phi_s[3]:

  detJ_s := simplify(detJ(3/4, 1/2)):
  
  for i from 0 to n_functions-1 do
    for j from 0 to n_functions-1 do
      # split integral over sector into the two parts
      mass_matrix[i,j] := mass_matrix[i,j] +
                          int(int(phi_s[i] * phi_s[j],
                                  y=0..x),
			      x=0..1/2) * detJ_s:
      mass_matrix[i,j] := mass_matrix[i,j] +
                          int(int(phi_s[i] * phi_s[j],
                                  y=0..1-x),
			      x=1/2..1) * detJ_s:
    od:
  od:  

  # sector 3
  phi_s[0] := 0:
  phi_s[1] := 0:
  phi_s[2] := x+y-1:
  phi_s[3] := y-x:
  phi_s[4] := 1 - phi_s[0] - phi_s[1] - phi_s[2] - phi_s[3]:

  detJ_s := simplify(detJ(1/2, 3/4)):
  
  for i from 0 to n_functions-1 do
    for j from 0 to n_functions-1 do
      # split integral over sector into the two parts
      mass_matrix[i,j] := mass_matrix[i,j] +
                          int(int(phi_s[i] * phi_s[j],
                                  y=0..x),
			      x=0..1/2) * detJ_s:
      mass_matrix[i,j] := mass_matrix[i,j] +
                          int(int(phi_s[i] * phi_s[j],
                                  y=0..1-x),
			      x=1/2..1) * detJ_s:
    od:
  od:  

  # sector 4
  phi_s[0] := 1-x-y:
  phi_s[1] := 0:
  phi_s[2] := 0:
  phi_s[3] := y-x:
  phi_s[4] := 1 - phi_s[0] - phi_s[1] - phi_s[2] - phi_s[3]:

  detJ_s := simplify(detJ(1/4, 1/2)):
  
  for i from 0 to n_functions-1 do
    for j from 0 to n_functions-1 do
      # split integral over sector into the two parts
      mass_matrix[i,j] := mass_matrix[i,j] +
                          int(int(phi_s[i] * phi_s[j],
                                  y=0..x),
			      x=0..1/2) * detJ_s:
      mass_matrix[i,j] := mass_matrix[i,j] +
                          int(int(phi_s[i] * phi_s[j],
                                  y=0..1-x),
			      x=1/2..1) * detJ_s:
    od:
  od:  
  
  print ("writing data to files"):
  readlib(C):
  C(prolongation, filename=prolongation_2d):
  C(restriction, filename=restriction_2d):
  C(array(1..2, [x[4], y[4]]), optimized, filename=crosspoint_2d):
  C(mass_matrix, optimized, filename=massmatrix_2d):
  
  --------------------------------------------------------------------
  
  Postprocess the files by the commands
  
  perl -pi -e 's/\[(\d+)\]\[(\d+)\]\[(\d+)\]/[$1]($2,$3)/g;' prolongation_2d
  perl -pi -e 's/.*= 0.0;\n//g;' prolongation_2d
  perl -pi -e 's/\[(\d+)\]\[(\d+)\]/($1,$2)/g;' massmatrix_2d
  perl -pi -e 's/(t\d+) =/const double $1 =/g;' massmatrix_2d

-----------------------------------------------------------------------------*/







#if deal_II_dimension == 1

template <>
FECrissCross<1>::FECrissCross () :
		FiniteElement<1> (0,0,0,0)
{
  Assert (false, ExcNotUseful());
};



template <>
double FECrissCross<1>::shape_value (const unsigned int, const Point<1> &) const {
  Assert (false, ExcNotUseful());
  return 0;
};



template <>
Tensor<1,1> FECrissCross<1>::shape_grad (const unsigned int, const Point<1> &) const {
  Assert (false, ExcNotUseful());
  return Point<1>();
};



template <>
Tensor<2,1> FECrissCross<1>::shape_grad_grad (const unsigned int, const Point<1> &) const {
  Assert (false, ExcNotUseful());
  return Tensor<2,1>();
};



template <>
void FECrissCross<1>::get_unit_support_points (vector<Point<1> >&) const {
  Assert (false, ExcNotUseful());
};



template <>
void FECrissCross<1>::get_support_points (const DoFHandler<1>::cell_iterator &,
					  const Boundary<1> &,
					  vector<Point<1> > &) const {
  Assert (false, ExcNotUseful());
};



template <>
void FECrissCross<1>::get_face_support_points (const DoFHandler<1>::face_iterator &,
					       const Boundary<1> &,
					       vector<Point<1> > &) const {
  Assert (false, ExcNotUseful());
};



template <>
void FECrissCross<1>::get_local_mass_matrix (const DoFHandler<1>::cell_iterator &,
					     const Boundary<1> &,
					     dFMatrix &) const {
  Assert (false, ExcNotUseful());
};



template <>
double  FECrissCross<1>::shape_value_transform (const unsigned int,
						const Point<1> &) const {
  Assert (false, ExcNotUseful());
  return 0;
};



template <>
Tensor<1,1> FECrissCross<1>::shape_grad_transform (const unsigned int,
						   const Point<1> &) const {
  Assert (false, ExcNotUseful());
  return Point<1>();
};



template <>
void FECrissCross<1>::get_face_jacobians (const DoFHandler<1>::face_iterator &,
					  const Boundary<1>       &,
					  const vector<Point<0> > &,
					  vector<double>      &) const {
  Assert (false, ExcNotUseful());
};



template <>
void FECrissCross<1>::get_subface_jacobians (const DoFHandler<1>::face_iterator &,
					     const unsigned int,
					     const vector<Point<0> > &,
					     vector<double>      &) const {
  Assert (false, ExcNotUseful());
};



template <>
void FECrissCross<1>::get_normal_vectors (const DoFHandler<1>::cell_iterator &,
					  const unsigned int,
					  const Boundary<1>       &,
					  const vector<Point<0> > &,
					  vector<Point<1> >       &) const {
  Assert (false, ExcNotUseful());
};



template <>
void FECrissCross<1>::get_normal_vectors (const DoFHandler<1>::cell_iterator &,
					  const unsigned int,
					  const unsigned int,
					  const vector<Point<0> > &,
					  vector<Point<1> >       &) const {
  Assert (false, ExcNotUseful());
};



template <>
void FECrissCross<1>::fill_fe_values (const DoFHandler<1>::cell_iterator &,
				      const vector<Point<1> >            &,
				      vector<dFMatrix>    &,
				      const bool,
				      vector<Point<1> > &,
				      const bool,
				      vector<Point<1> > &,
				      const bool,
				      const dFMatrix      &,
				      const vector<vector<Tensor<1,1> > > &,
				      const Boundary<1> &) const {
  Assert (false, ExcNotUseful());
};

#endif






#if deal_II_dimension == 2

template <>
FECrissCross<2>::FECrissCross () :
		FiniteElement<2> (1,0,1,5)
{
  interface_constraints(0,0) = 1./2.;
  interface_constraints(0,1) = 1./2.;

  prolongation[0](0,0) = 1.0;
  prolongation[0](1,0) = 1.0/2.0;
  prolongation[0](1,1) = 1.0/2.0;
  prolongation[0](2,4) = 1.0;
  prolongation[0](3,0) = 1.0/2.0;
  prolongation[0](3,3) = 1.0/2.0;
  prolongation[0](4,0) = 1.0/2.0;
  prolongation[0](4,4) = 1.0/2.0;
  prolongation[1](0,0) = 1.0/2.0;
  prolongation[1](0,1) = 1.0/2.0;
  prolongation[1](1,1) = 1.0;
  prolongation[1](2,1) = 1.0/2.0;
  prolongation[1](2,2) = 1.0/2.0;
  prolongation[1](3,4) = 1.0;
  prolongation[1](4,1) = 1.0/2.0;
  prolongation[1](4,4) = 1.0/2.0;
  prolongation[2](0,4) = 1.0;
  prolongation[2](1,1) = 1.0/2.0;
  prolongation[2](1,2) = 1.0/2.0;
  prolongation[2](2,2) = 1.0;
  prolongation[2](3,2) = 1.0/2.0;
  prolongation[2](3,3) = 1.0/2.0;
  prolongation[2](4,2) = 1.0/2.0;
  prolongation[2](4,4) = 1.0/2.0;
  prolongation[3](0,0) = 1.0/2.0;
  prolongation[3](0,3) = 1.0/2.0;
  prolongation[3](1,4) = 1.0;
  prolongation[3](2,2) = 1.0/2.0;
  prolongation[3](2,3) = 1.0/2.0;
  prolongation[3](3,3) = 1.0;
  prolongation[3](4,3) = 1.0/2.0;
  prolongation[3](4,4) = 1.0/2.0;

  restriction[0](0,0) = 1.0;
  restriction[0](4,2) = 1.0;
  restriction[1](1,1) = 1.0;
  restriction[1](4,3) = 1.0;
  restriction[2](2,2) = 1.0;
  restriction[2](4,0) = 1.0;
  restriction[3](3,3) = 1.0;
  restriction[3](4,1) = 1.0;
};



template <>
inline
double FECrissCross<2>::shape_value (const unsigned int i,
				     const Point<2>    &p) const {
  Assert((i<total_dofs), ExcInvalidIndex(i));

  const double x = p(0),
	       y = p(1);
  switch (i)
    {
      case 0: return ((y<1-x) ? 1-x-y : 0);
      case 1: return ((y<x)   ? x-y   : 0);
      case 2: return ((y>1-x) ? x+y-1 : 0);
      case 3: return ((y>x)   ? y-x   : 0);

					     // I am too lazy to optimize the
					     // following myself. Let the
					     // compiler do this
      case 4: return (1-(((y<1-x) ? 1-x-y : 0) +
			 ((y<x)   ? x-y   : 0) +
			 ((y>1-x) ? x+y-1 : 0) +
			 ((y>x)   ? y-x   : 0)));
    }
  return 0.;
};



template <>
inline
Tensor<1,2> FECrissCross<2>::shape_grad (const unsigned int i, const Point<2> &p) const {
  Assert((i<total_dofs), ExcInvalidIndex(i));

  const double x = p(0),
	       y = p(1);  
  switch (i)
    {
      case 0: return ((y<1-x) ? Point<2>(-1,-1) : Point<2>(0,0));
      case 1: return ((y<x)   ? Point<2>(1,-1)  : Point<2>(0,0));
      case 2: return ((y>1-x) ? Point<2>(1,1)   : Point<2>(0,0));
      case 3: return ((y>x)   ? Point<2>(-1,1)  : Point<2>(0,0));

      					     // I am too lazy to optimize the
					     // following myself. Let the
					     // compiler do this
      case 4: return -1.*(((y<1-x) ? Point<2>(-1,-1) : Point<2>(0,0)) +
			  ((y<x)   ? Point<2>(1,-1)  : Point<2>(0,0)) +
			  ((y>1-x) ? Point<2>(1,1)   : Point<2>(0,0)) +
			  ((y>x)   ? Point<2>(-1,1)  : Point<2>(0,0)));
    }
  return Point<2>();
};


template <>
inline
Tensor<2,2>
FECrissCross<2>::shape_grad_grad (const unsigned int i,
				  const Point<2> &) const
{
  Assert((i<total_dofs), ExcInvalidIndex(i));
				   // second derivatives on the unit cell
				   // are always zero, at least almost
				   // everywhere. see the doc for more
				   // info
  return Tensor<2,2>();
};



template <>
void FECrissCross<2>::get_unit_support_points (vector<Point<2> > &unit_points) const {
  Assert(unit_points.size()==total_dofs,
	  ExcWrongFieldDimension (unit_points.size(), total_dofs));

  unit_points[0] = Point<2> (0,0);
  unit_points[1] = Point<2> (1,0);
  unit_points[2] = Point<2> (1,1);
  unit_points[3] = Point<2> (0,1);
  unit_points[4] = Point<2> (0.5, 0.5);
};



template <>
void FECrissCross<2>::get_support_points (const DoFHandler<2>::cell_iterator &cell,
					  const Boundary<2> &,
					  vector<Point<2> > &support_points) const {
  const unsigned int dim = 2;
  
  Assert (support_points.size() == total_dofs,
	  ExcWrongFieldDimension (support_points.size(), total_dofs));

				   // copy vertices
  for (unsigned int vertex=0; vertex<GeometryInfo<dim>::vertices_per_cell; ++vertex)
    support_points[vertex] = cell->vertex(vertex);

/*
  last support point is the common point of the two diagonals; the formula for
  the computation is a bit lengthy but straightforward. You can get it with
  the small Maple script printed at the beginning of this file.
*/
  const double x0 = cell->vertex(0)(0),
	       y0 = cell->vertex(0)(1),
	       x1 = cell->vertex(1)(0),
	       y1 = cell->vertex(1)(1),
	       x2 = cell->vertex(2)(0),
	       y2 = cell->vertex(2)(1),
	       x3 = cell->vertex(3)(0),
	       y3 = cell->vertex(3)(1);
  const double t1 = x0*y1;
  const double t2 = x0*y3;
  const double t4 = x1*y0;
  const double t5 = x3*y0;
  const double t14 = (t1-t2+x1*y3-t4+t5-x3*y1)/(t1-t2-x2*y1+x2*y3-t4+x1*y2+t5-x3*y2);
  const double t15 = 1.0-t14;
  support_points[4](0) = t15*x0+t14*x2;
  support_points[4](1) = t15*y0+t14*y2;
};



template <>
void FECrissCross<2>::get_face_support_points (const DoFHandler<2>::face_iterator &face,
					       const Boundary<2> &,
					       vector<Point<2> > &support_points) const {
  const unsigned int dim = 2;
  
  Assert ((support_points.size() == dofs_per_face) &&
	  (support_points.size() == GeometryInfo<dim>::vertices_per_face),
	  ExcWrongFieldDimension (support_points.size(),
				  GeometryInfo<dim>::vertices_per_face));

  for (unsigned int vertex=0; vertex<dofs_per_face; ++vertex)
    support_points[vertex] = face->vertex(vertex);
};



template <>
void FECrissCross<2>::get_local_mass_matrix (const DoFHandler<2>::cell_iterator &cell,
					     const Boundary<2> &,
					     dFMatrix &mass_matrix) const {
  Assert (mass_matrix.n() == total_dofs,
	  ExcWrongFieldDimension(mass_matrix.n(),total_dofs));
  Assert (mass_matrix.m() == total_dofs,
	  ExcWrongFieldDimension(mass_matrix.m(),total_dofs));

  const double x[4] = { cell->vertex(0)(0),
			cell->vertex(1)(0),
			cell->vertex(2)(0),
			cell->vertex(3)(0)  };
  const double y[4] = { cell->vertex(0)(1),
			cell->vertex(1)(1),
			cell->vertex(2)(1),
			cell->vertex(3)(1)  };

  const double t1 = x[3]*x[2];
  const double t2 = y[1]*y[1];
  const double t5 = x[0]*x[0];
  const double t7 = t5*y[3]*y[2];
  const double t8 = x[3]*x[1];
  const double t9 = y[0]*y[0];
  const double t10 = t8*t9;
  const double t11 = t5*y[1];
  const double t12 = t11*y[3];
  const double t13 = t11*y[2];
  const double t14 = x[1]*x[1];
  const double t15 = t14*y[0];
  const double t19 = x[2]*x[1]*t9;
  const double t20 = t1*t9;
  const double t24 = x[0]*x[2];
  const double t25 = t24*t2;
  const double t29 = t15*y[2];
  const double t30 = x[3]*y[0];
  const double t31 = x[1]*y[2];
  const double t32 = t30*t31;
  const double t33 = x[2]*y[3];
  const double t34 = x[1]*y[0];
  const double t35 = t33*t34;
  const double t37 = x[2]*y[1];
  const double t38 = t37*t34;
  const double t39 = y[1]*x[1];
  const double t42 = -2.0*t1*t2-t7+t10+t12+t13+2.0*t15*y[3]+t19-t20-2.0*t14*y[2]*y[3]+
		     t25+2.0*x[0]*t2*x[3]+t29-t32-2.0*t35-t38-2.0*t39*t30;
  const double t43 = y[3]*x[1];
  const double t46 = x[3]*y[2];
  const double t49 = t37*t30;
  const double t51 = x[0]*y[1];
  const double t52 = t51*t46;
  const double t54 = t51*t33;
  const double t55 = x[0]*y[0];
  const double t56 = t55*t37;
  const double t59 = t51*t30;
  const double t60 = t51*t31;
  const double t61 = t55*t31;
  const double t62 = t55*t46;
  const double t63 = x[0]*y[3];
  const double t64 = t63*t34;
  const double t65 = t55*t33;
  const double t66 = t63*t31;
  const double t72 = 2.0*t37*t43+2.0*t39*t46+3.0*t49-2.0*t52-t54-t56-2.0*t51*t43-t59-t60
		     -t61+t62-t64+t65+3.0*t66-t14*t9-t5*t2+2.0*t51*t34;
  const double t75 = 1/(t51-t63-t37+t33-t34+t31+t30-t46);
  const double t76 = (t42+t72)*t75;
  const double t77 = y[3]*y[3];
  const double t81 = x[3]*x[3];
  const double t82 = y[1]*t81;
  const double t86 = t81*y[0]*y[2];
  const double t90 = t24*t77;
  const double t94 = -t7-t10-t12-2.0*x[0]*t77*x[1]+t13+t19+2.0*t82*y[2]-t86-t20+t81*t9
		     -2.0*t82*y[0]-t90+t32-3.0*t35+2.0*t49-3.0*t52;
  const double t96 = t63*t46;
  const double t97 = t30*t33;
  const double t98 = x[3]*y[3];
  const double t114 = t54-t56+t59-t61+t62+t64+t65+2.0*t66+t96+t97+2.0*t51*t98-2.0*t63*
		      t30-2.0*y[1]*x[3]*t33+2.0*t98*t34-2.0*t98*t31+2.0*t77*x[1]*x[2]+t5*t77;
  const double t116 = (t94+t114)*t75;
  const double t118 = t76/24.0;
  const double t119 = t116/24.0;
  const double t121 = -t118+t116/8.0;
  const double t122 = x[0]*y[2];
  const double t123 = t122*t37;
  const double t124 = t122*t33;
  const double t125 = y[2]*y[2];
  const double t126 = x[0]*t125;
  const double t127 = t126*x[1];
  const double t128 = t126*x[3];
  const double t129 = x[2]*x[2];
  const double t131 = y[1]*t129;
  const double t132 = t131*y[0];
  const double t133 = t131*y[3];
  const double t134 = -t25+t60+t123+t54-t52-t124-t127+t128+t129*t2-t132-t133;
  const double t135 = t37*t46;
  const double t139 = x[2]*y[0];
  const double t140 = t139*t46;
  const double t141 = t8*t125;
  const double t143 = t129*y[0]*y[3];
  const double t144 = t33*t31;
  const double t145 = t139*t31;
  const double t146 = t38+t135-2.0*t37*t31-t29+t14*t125-t140-t35-t141+t143+t144+t32+t145
		      ;
  const double t148 = (t134+t146)*t75;
  const double t150 = t148/24.0;
  const double t152 = -t118+t148/8.0;
  const double t153 = t123-t54+t66+t128-t96-t124-t127+t90+t49-t132-t135;
  const double t158 = t133+2.0*t33*t46+t86-t81*t125-t129*t77+t141-t97-t32-t144-t140+t143
		      +t145;
  const double t160 = (t153+t158)*t75;
  const double t162 = t160/24.0;
  const double t164 = 7.0/24.0*t160;
  const double t165 = -5.0/24.0*t148+t164;
  const double t168 = t164-5.0/24.0*t116;
  mass_matrix(0,0) = -t76/12.0+t116/12.0;
  mass_matrix(0,1) = -t118;
  mass_matrix(0,2) = 0.0;
  mass_matrix(0,3) = -t119;
  mass_matrix(0,4) = t121;
  mass_matrix(1,0) = -t118;
  mass_matrix(1,1) = -t76/12.0+t148/12.0;
  mass_matrix(1,2) = -t150;
  mass_matrix(1,3) = 0.0;
  mass_matrix(1,4) = t152;
  mass_matrix(2,0) = 0.0;
  mass_matrix(2,1) = -t150;
  mass_matrix(2,2) = t148/12.0-t160/12.0;
  mass_matrix(2,3) = -t162;
  mass_matrix(2,4) = t165;
  mass_matrix(3,0) = -t119;
  mass_matrix(3,1) = 0.0;
  mass_matrix(3,2) = -t162;
  mass_matrix(3,3) = -t160/12.0+t116/12.0;
  mass_matrix(3,4) = t168;
  mass_matrix(4,0) = t121;
  mass_matrix(4,1) = t152;
  mass_matrix(4,2) = t165;
  mass_matrix(4,3) = t168;
  mass_matrix(4,4) = -t76/12.0+7.0/12.0*t148-17.0/12.0*t160+7.0/12.0*t116;
};



template <>
inline
double FECrissCross<2>::shape_value_transform (const unsigned int i,
					       const Point<2> &p) const {
				   // use an isoparametric ansatz
  return shape_value(i,p);
};



template <>
Tensor<1,2> FECrissCross<2>::shape_grad_transform (const unsigned int i,
						   const Point<2> &p) const {
				   // use an isoparametric ansatz
  return shape_grad(i,p);  
};



template <>
void FECrissCross<2>::get_face_jacobians (const DoFHandler<2>::face_iterator &face,
					  const Boundary<2>       &,
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
void FECrissCross<2>::get_subface_jacobians (const DoFHandler<2>::face_iterator &face,
					     const unsigned int,
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
void FECrissCross<2>::get_normal_vectors (const DoFHandler<2>::cell_iterator &cell,
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
void FECrissCross<2>::get_normal_vectors (const DoFHandler<2>::cell_iterator &cell,
					  const unsigned int       face_no,
					  const unsigned int,
					  const vector<Point<1> > &unit_points,
					  vector<Point<2> >       &normal_vectors) const {
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



template <int dim>
void FECrissCross<dim>::fill_fe_values (const DoFHandler<dim>::cell_iterator &cell,
					const vector<Point<dim> >            &unit_points,
					vector<Tensor<2,dim> >               &jacobians,
					const bool              compute_jacobians,
					vector<Tensor<3,dim> > &jacobians_grad,
					const bool              compute_jacobians_grad,
					vector<Point<dim> >    &support_points,
					const bool,
					vector<Point<dim> >    &q_points,
					const bool              compute_q_points,
					const dFMatrix         &shape_values_transform,
					const vector<vector<Tensor<1,dim> > > &/*shape_grad_transform*/,
					const Boundary<dim> &boundary) const {
  Assert (jacobians.size() == unit_points.size(),
	  ExcWrongFieldDimension(jacobians.size(), unit_points.size()));
  Assert (q_points.size() == unit_points.size(),
	  ExcWrongFieldDimension(q_points.size(), unit_points.size()));
  Assert (support_points.size() == total_dofs,
	  ExcWrongFieldDimension(support_points.size(), total_dofs));

  
  unsigned int n_points=unit_points.size();

				   // we need the support points in any
				   // way, wanted or not by the user
  get_support_points (cell, boundary, support_points);

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
      for (unsigned int j=0; j<n_transform_functions; ++j) 
	for (unsigned int l=0; l<n_points; ++l) 
	  q_points[l] += support_points[j] * shape_values_transform(j, l);
    };
  

/* jacobi matrices: compute d(x)/d(xi) and invert this
   Let M(l) be the inverse of J at the quadrature point l, then
     M_{ij}(l) = sum_s p_i(s) d(N_s(l))/d(xi_j)
   where p_i(s) is the i-th coordinate of the s-th vertex vector,
   N_s(l) is the value of the s-th shape function at the
   quadrature point l.

   We could therefore write:
   l=0..n_points-1
     i=0..dim-1
       j=0..dim-1
         M_{ij}(l) = 0
	 s=0..n_support_points
	   M_{ij}(l) += p_i(s) d(N_s(l))/d(xi_j)

  However, we rewrite the loops to only compute the gradient once for
  each integration point and basis function.

  The scheme laid down above was originally used. Due to recent advances
  in the authors understanding of most basic things, it was dropped and
  replaced by the following version. See #FELinearMapping<dim>::fill_fe_values#
  for more information on this.
*/

  Point<dim> vertices[GeometryInfo<dim>::vertices_per_cell];
  for (unsigned int l=0; l<GeometryInfo<dim>::vertices_per_cell; ++l)
    vertices[l] = cell->vertex(l);

  if (compute_jacobians)
    switch (dim)
      {
	case 1:
	      for (unsigned int point=0; point<n_points; ++point)
		jacobians[point][0][0] = 1./(vertices[1](0)-vertices[0](0));
	      break;
		
	case 2:
	{
	  for (unsigned int point=0; point<n_points; ++point)
	    {	    
	      const double xi = unit_points[point](0);
	      const double eta= unit_points[point](1);
	
	      const double t6 = vertices[0](0)*vertices[3](1);
	      const double t8 = vertices[2](0)*xi;
	      const double t10 = vertices[1](0)*eta;
	      const double t12 = vertices[3](0)*vertices[1](1);
	      const double t16 = vertices[3](0)*xi;
	      const double t20 = vertices[0](0)*vertices[1](1);
	      const double t22 = vertices[0](0)*vertices[2](1);
	      const double t24 = t6*xi-t8*vertices[1](1)-t10*vertices[3](1)+
				 t12*eta-vertices[3](0)*vertices[2](1)*eta-
				 t16*vertices[0](1)+t16*vertices[1](1)-t12+
				 vertices[3](0)*vertices[0](1)-t20*eta+t22*eta;
	      const double t28 = vertices[1](0)*vertices[3](1);
	      const double t31 = vertices[2](0)*eta;
	      const double t36 = t8*vertices[0](1)+vertices[1](0)*vertices[2](1)*xi-
				 t28*xi+t10*vertices[0](1)-t31*vertices[0](1)+
				 t31*vertices[3](1)+t20-t6-vertices[1](0)*
				 vertices[0](1)+t28-t22*xi;
	      const double t38 = 1/(t24+t36);

	      jacobians[point][0][0] = (-vertices[0](1)+vertices[0](1)*xi-
					vertices[1](1)*xi+vertices[2](1)*xi+
					vertices[3](1)-vertices[3](1)*xi)*t38;
	      jacobians[point][0][1] = -(-vertices[0](0)+vertices[0](0)*xi-
					 vertices[1](0)*xi+t8+vertices[3](0)-t16)*t38;
	      jacobians[point][1][0] = -(-vertices[0](1)+vertices[0](1)*eta+
					 vertices[1](1)-vertices[1](1)*eta+
					 vertices[2](1)*eta-vertices[3](1)*eta)*t38;
	      jacobians[point][1][1] = (-vertices[0](0)+vertices[0](0)*eta+
					vertices[1](0)-t10+t31-vertices[3](0)*eta)*t38;
	    };
	  
	  break;
	};

	default:
					       // not implemented at present
	      Assert (false, ExcNotImplemented());
      };


  if (compute_jacobians_grad)
    switch (dim) 
      {
	case 1:
	{
					   // derivative of the
					   // jacobian is always zero
					   // for a linear mapping in 1d
	  for (unsigned int point=0; point<n_points; ++point)
	    jacobians_grad[point][0][0][0] = 0;
	  break;
	};

	case 2:
	{
	  for (unsigned int point=0; point<n_points; ++point)
	    {
	      const double xi = unit_points[point](0);
	      const double eta= unit_points[point](1);
	
	      const double t2 = vertices[1](0)*eta;
	      const double t4 = vertices[3](0)*vertices[2](1);
	      const double t6 = vertices[0](0)*vertices[2](1);
	      const double t8 = vertices[0](0)*vertices[3](1);
	      const double t10 = vertices[3](0)*xi;
	      const double t13 = vertices[2](0)*xi;
	      const double t16 = vertices[3](0)*vertices[1](1);
	      const double t18 = vertices[0](0)*vertices[1](1);
	      const double t19 = vertices[3](0)*vertices[0](1);
	      const double t20 = -t2*vertices[3](1)-t4*eta-t6*xi+t8*xi-
				 t10*vertices[0](1)+t10*vertices[1](1)+
				 t13*vertices[0](1)-t13*vertices[1](1)+t16
				 *eta+t18+t19;
	      const double t23 = vertices[1](0)*vertices[3](1);
	      const double t26 = vertices[2](0)*eta;
	      const double t29 = vertices[1](0)*vertices[0](1);
	      const double t30 = vertices[1](0)*vertices[2](1);
	      const double t32 = -t16-t18*eta+t6*eta-t23*xi+t2*vertices[0](1)-
				 t26*vertices[0](1)+t26*vertices[3](1)-
				 t8-t29+t23+t30
				 *xi;
	      const double t33 = t20+t32;
	      const double t34 = 1/t33;
	      const double t35 = (vertices[0](1)-vertices[1](1)+
				  vertices[2](1)-vertices[3](1))*t34;
	      const double t41 = t33*t33;
	      const double t42 = 1/t41;
	      const double t43 = (-vertices[0](1)+vertices[0](1)*xi-
				  vertices[1](1)*xi+vertices[2](1)*xi+
				  vertices[3](1)-vertices[3](1)*xi)*t42;
	      const double t44 = vertices[2](0)*vertices[0](1);
	      const double t46 = -t6+t8-t19+t16+t44-
				 vertices[2](0)*vertices[1](1)-
				 t23+t30;
	      const double t50 = (vertices[0](0)-vertices[1](0)+
				  vertices[2](0)-vertices[3](0))*t34;
	      const double t54 = (-vertices[0](0)+vertices[0](0)*xi-
				  vertices[1](0)*xi+t13+
				  vertices[3](0)-t10)*t42;
	      const double t62 = (-vertices[0](1)+vertices[0](1)*eta+
				  vertices[1](1)-vertices[1](1)*eta+
				  vertices[2](1)*eta-
				  vertices[3](1)*eta)*t42;
	      const double t67 = (-vertices[0](0)+vertices[0](0)*eta+
				  vertices[1](0)-t2+t26-vertices[3](0)*eta)*t42;
	      const double t70 = -t23-t4+t16-t18+t6+t29-t44+
				 vertices[2](0)*vertices[3](1);
	      jacobians_grad[point][0][0][0] = t35-t43*t46;
	      jacobians_grad[point][0][0][1] = -t50+t54*t46;
	      jacobians_grad[point][0][1][0] = t62*t46;
	      jacobians_grad[point][0][1][1] = -t67*t46;
	      jacobians_grad[point][1][0][0] = -t43*t70;
	      jacobians_grad[point][1][0][1] = t54*t70;
	      jacobians_grad[point][1][1][0] = -t35+t62*t70;
	      jacobians_grad[point][1][1][1] = t50-t67*t70;
	    };
	  break;
	  
	};
	
	default:
					       // not implemented at present
	      Assert (false, ExcNotImplemented());
      };	      
};

#endif





/*--------------------------- QCrissCross* ------------------------------------*/


#if deal_II_dimension == 1

template <>
QCrissCross1<1>::QCrissCross1 () :
		Quadrature<1> (1)
{
  Assert (false, ExcNotUseful());
};

#endif



#if deal_II_dimension == 2

template <>
QCrissCross1<2>::QCrissCross1 () :
		Quadrature<2> (4)
{
				   // let quadrature points be the
				   // barycenters of the four triangles
  quadrature_points[0] = Point<2>(1./2., 1./6.);
  quadrature_points[1] = Point<2>(5./6., 1./2.);
  quadrature_points[2] = Point<2>(1./2., 5./6.);
  quadrature_points[3] = Point<2>(1./6., 1./2.);

  weights[0] = 1./4.;
  weights[1] = 1./4.;
  weights[2] = 1./4.;
  weights[3] = 1./4.;
};

#endif



// explicit instantiations

template class FECrissCross<deal_II_dimension>;
template class QCrissCross1<deal_II_dimension>;
