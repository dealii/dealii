/* $Id$ */

#include <fe/fe_lib.criss_cross.h>
#include <grid/tria_iterator.h>
#include <grid/dof_accessor.h>
#include <grid/geometry_info.h>




/*-----------------------------------2d------------------------------------
  Maple script to automate some of the error-prone computations on
  this strange sort of element.

  n_functions      := 5:

  # note: ansatz_points[i] is a vector which is indexed from
  # one and not from zero!
  ansatz_points[0] := [0,0]:
  ansatz_points[1] := [1,0]:
  ansatz_points[2] := [1,1]:
  ansatz_points[3] := [0,1]:
  ansatz_points[4] := [1/2,1/2]:
  
  phi[0] := proc(x,y) if(y<1-x) then 1-x-y; else 0; fi; end:
  phi[1] := proc(x,y) if(y<x)   then x-y;   else 0; fi; end:
  phi[2] := proc(x,y) if(y>1-x) then x+y-1; else 0; fi; end:
  phi[3] := proc(x,y) if(y>x)   then y-x;   else 0; fi; end:
  phi[4] := proc(x,y) 1 - phi[0](x,y) - phi[1](x,y)
                        - phi[2](x,y) - phi[3](x,y) ; end:

  #points on children: let them be indexed one-based, as are
  #the ansatz_points
  points[0] := array(0..n_functions-1, 1..2):
  points[1] := array(0..n_functions-1, 1..2):
  points[2] := array(0..n_functions-1, 1..2):
  points[3] := array(0..n_functions-1, 1..2):
  for i from 0 to n_functions-1 do
    points[0][i,1] := ansatz_points[i][1]/2:
    points[0][i,2] := ansatz_points[i][2]/2:
    
    points[1][i,1] := ansatz_points[i][1]/2+1/2:
    points[1][i,2] := ansatz_points[i][2]/2:

    points[2][i,1] := ansatz_points[i][1]/2+1/2:
    points[2][i,2] := ansatz_points[i][2]/2+1/2:

    points[3][i,1] := ansatz_points[i][1]/2:
    points[3][i,2] := ansatz_points[i][2]/2+1/2:
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

  eq_sys := {(1-t)*x0 + t*x2 = (1-s)*x1 + s*x3,
             (1-t)*y0 + t*y2 = (1-s)*y1 + s*y3}:
  solution := solve (eq_sys, {s,t});

  xs := subs (solution, (1-t)*x0 + t*x2):
  ys := subs (solution, (1-t)*y0 + t*y2):
  ps := array(1..2, [xs, ys]):

  print ("writing data to files"):
  readlib(C):
  C(prolongation, filename=prolongation_2d):
  C(ps, filename=crosspoint):

  --------------------------------------------------------------------
  
  Postprocess the prolongation matrix by the commands
  
  perl -pi -e 's/\[(\d+)\]\[(\d+)\]\[(\d+)\]/[$1]($2,$3)/g;' prolongation_2d
  perl -pi -e 's/.*= 0.0;\n//g;' prolongation_2d

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
Point<1> FECrissCross<1>::shape_grad (const unsigned int, const Point<1> &) const {
  Assert (false, ExcNotUseful());
  return Point<1>();
};



template <>
void FECrissCross<1>::get_unit_ansatz_points (vector<Point<1> >&) const {
  Assert (false, ExcNotUseful());
};



template <>
void FECrissCross<1>::get_ansatz_points (const DoFHandler<1>::cell_iterator &,
					 const Boundary<1> &,
					 vector<Point<1> > &) const {
  Assert (false, ExcNotUseful());
};



template <>
void FECrissCross<1>::get_face_ansatz_points (const DoFHandler<1>::face_iterator &,
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
Point<1> FECrissCross<1>::shape_grad_transform (const unsigned int,
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
				      const vector<vector<Point<1> > > &,
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
Point<2> FECrissCross<2>::shape_grad (const unsigned int i, const Point<2> &p) const {
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
void FECrissCross<2>::get_unit_ansatz_points (vector<Point<2> > &unit_points) const {
  Assert(unit_points.size()==total_dofs,
	  ExcWrongFieldDimension (unit_points.size(), total_dofs));

  unit_points[0] = Point<2> (0,0);
  unit_points[1] = Point<2> (1,0);
  unit_points[2] = Point<2> (1,1);
  unit_points[3] = Point<2> (0,1);
  unit_points[4] = Point<2> (0.5, 0.5);
};



template <>
void FECrissCross<2>::get_ansatz_points (const DoFHandler<2>::cell_iterator &cell,
					 const Boundary<2> &,
					 vector<Point<2> > &ansatz_points) const {
  const unsigned int dim = 2;
  
  Assert (ansatz_points.size() == total_dofs,
	  ExcWrongFieldDimension (ansatz_points.size(), total_dofs));

				   // copy vertices
  for (unsigned int vertex=0; vertex<GeometryInfo<dim>::vertices_per_cell; ++vertex)
    ansatz_points[vertex] = cell->vertex(vertex);

/*
  last ansatz point is the common point of the two diagonals; the formula for
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
  ansatz_points[4](0) = t15*x0+t14*x2;
  ansatz_points[4](1) = t15*y0+t14*y2;
};



template <>
void FECrissCross<2>::get_face_ansatz_points (const DoFHandler<2>::face_iterator &face,
					      const Boundary<2> &,
					      vector<Point<2> > &ansatz_points) const {
  const unsigned int dim = 2;
  
  Assert ((ansatz_points.size() == dofs_per_face) &&
	  (ansatz_points.size() == GeometryInfo<dim>::vertices_per_face),
	  ExcWrongFieldDimension (ansatz_points.size(),
				  GeometryInfo<dim>::vertices_per_face));

  for (unsigned int vertex=0; vertex<dofs_per_face; ++vertex)
    ansatz_points[vertex] = face->vertex(vertex);
};



template <>
void FECrissCross<2>::get_local_mass_matrix (const DoFHandler<2>::cell_iterator &,
					     const Boundary<2> &,
					     dFMatrix &local_mass_matrix) const {
  Assert (local_mass_matrix.n() == total_dofs,
	  ExcWrongFieldDimension(local_mass_matrix.n(),total_dofs));
  Assert (local_mass_matrix.m() == total_dofs,
	  ExcWrongFieldDimension(local_mass_matrix.m(),total_dofs));

				   // in this special element, some of the
				   // entries are zero (which is not the
				   // case for most other elements, so
				   // we first reset all elements and only
				   // fill in those that are nonzero
  local_mass_matrix.clear ();
  
  Assert (false, ExcNotUseful());
};



template <>
inline
double FECrissCross<2>::shape_value_transform (const unsigned int i,
						const Point<2> &p) const {
				   // use an isoparametric ansatz
  return shape_value(i,p);
};



template <>
Point<2> FECrissCross<2>::shape_grad_transform (const unsigned int i,
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
					vector<dFMatrix>    &jacobians,
					const bool           compute_jacobians,
					vector<Point<dim> > &ansatz_points,
					const bool,
					vector<Point<dim> > &q_points,
					const bool           compute_q_points,
					const dFMatrix      &shape_values_transform,
					const vector<vector<Point<dim> > > &shape_grad_transform,
					const Boundary<dim> &boundary) const {
  Assert (jacobians.size() == unit_points.size(),
	  ExcWrongFieldDimension(jacobians.size(), unit_points.size()));
  Assert (q_points.size() == unit_points.size(),
	  ExcWrongFieldDimension(q_points.size(), unit_points.size()));
  Assert (ansatz_points.size() == total_dofs,
	  ExcWrongFieldDimension(ansatz_points.size(), total_dofs));

  
  unsigned int n_points=unit_points.size();

				   // we need the ansatz points in any
				   // way, wanted or not by the user
  get_ansatz_points (cell, boundary, ansatz_points);

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
	  q_points[l] += ansatz_points[j] * shape_values_transform(j, l);
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
	 s=0..n_ansatz_points
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
	  for (unsigned int s=0; s<n_transform_functions; ++s)
	    {
					       // we want the linear transform,
					       // so use that function
	      const Point<dim> gradient = shape_grad_transform[s][l];
	      for (unsigned int i=0; i<dim; ++i)
		for (unsigned int j=0; j<dim; ++j)
		  M(i,j) += ansatz_points[s](i) * gradient(j);
	    };
	  jacobians[l].invert(M);
	};
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
