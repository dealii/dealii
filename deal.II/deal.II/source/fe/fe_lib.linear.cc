/* $Id$ */

#include <fe/fe_lib.h>
#include <grid/tria_iterator.h>
#include <grid/dof_accessor.h>
#include <grid/geometry_info.h>
#include <algorithm>




#if deal_II_dimension == 1

template <>
FELinear<1>::FELinear () :
		FiniteElement<1> (1, 0)
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
  restriction[0](0,1) = 1./2.;
  restriction[0](1,1) = 1./2.;

  restriction[1](0,0) = 1./2.;
  restriction[1](1,0) = 1./2.;
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
FELinear<1>::shape_value(const unsigned int i,
			 const Point<1>     &p) const
{
  Assert((i<total_dofs), ExcInvalidIndex(i));
  switch (i)
    {
    case 0: return 1.-p(0);
    case 1: return p(0);
    }
  return 0.;
}



template <>
inline
Point<1>
FELinear<1>::shape_grad(const unsigned int i,
			const Point<1>&) const
{
  Assert((i<total_dofs), ExcInvalidIndex(i));
  switch (i)
    {
    case 0: return Point<1>(-1.);
    case 1: return Point<1>(1.);
    }
  return Point<1>();
};



template <>
void FELinear<1>::fill_fe_values (const DoFHandler<1>::cell_iterator &cell,
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
void FELinear<1>::get_ansatz_points (const typename DoFHandler<1>::cell_iterator &cell,
				     const Boundary<1>  &boundary,
				     vector<Point<1> >  &ansatz_points) const {
  FiniteElement<1>::get_ansatz_points (cell, boundary, ansatz_points);
};



template <>
void FELinear<1>::get_face_ansatz_points (const typename DoFHandler<1>::face_iterator &,
					  const Boundary<1>  &,
					  vector<Point<1> >  &) const {
  Assert (false, ExcInternalError());
};



template <>
void FELinear<1>::get_face_jacobians (const DoFHandler<1>::face_iterator &,
				      const Boundary<1>         &,
				      const vector<Point<0> > &,
				      vector<double>      &) const {
  Assert (false, ExcInternalError());
};



template <>
void FELinear<1>::get_subface_jacobians (const DoFHandler<1>::face_iterator &,
					 const unsigned int           ,
					 const vector<Point<0> > &,
					 vector<double>      &) const {
  Assert (false, ExcInternalError());
};



template <>
void FELinear<1>::get_normal_vectors (const DoFHandler<1>::cell_iterator &,
				      const unsigned int,
				      const Boundary<1> &,
				      const vector<Point<0> > &,
				      vector<Point<1> > &) const {
  Assert (false, ExcInternalError());
};



template <>
void FELinear<1>::get_normal_vectors (const DoFHandler<1>::cell_iterator &,
				      const unsigned int,
				      const unsigned int,
				      const vector<Point<0> > &,
				      vector<Point<1> > &) const {
  Assert (false, ExcInternalError());
};



template <>
void FELinear<1>::get_local_mass_matrix (const DoFHandler<1>::cell_iterator &cell,
					 const Boundary<1> &,
					 dFMatrix &local_mass_matrix) const {
  Assert (local_mass_matrix.n() == total_dofs,
	  ExcWrongFieldDimension(local_mass_matrix.n(),total_dofs));
  Assert (local_mass_matrix.m() == total_dofs,
	  ExcWrongFieldDimension(local_mass_matrix.m(),total_dofs));

  const double h = cell->vertex(1)(0) - cell->vertex(0)(0);
  Assert (h>0, ExcJacobiDeterminantHasWrongSign());

  local_mass_matrix(0,0) = local_mass_matrix(1,1) = 1./3.*h;
  local_mass_matrix(0,1) = local_mass_matrix(1,0) = 1./6.*h;
};

#endif




#if deal_II_dimension == 2

template <>
FELinear<2>::FELinear () :
		FiniteElement<2> (1, 0, 0)
{
  interface_constraints(0,0) = 1./2.;
  interface_constraints(0,1) = 1./2.;

  restriction[0](0,0) = 1.0;
  restriction[0](0,1) = 1./2.;
  restriction[0](1,1) = 1./2.;
  restriction[0](0,3) = 1./2.;
  restriction[0](3,3) = 1./2.;
  restriction[0](0,2) = 1./4.;
  restriction[0](1,2) = 1./4.;
  restriction[0](2,2) = 1./4.;
  restriction[0](3,2) = 1./4.;

  restriction[1](1,1) = 1.0;
  restriction[1](0,0) = 1./2.;
  restriction[1](1,0) = 1./2.;
  restriction[1](1,2) = 1./2.;
  restriction[1](2,2) = 1./2.;
  restriction[1](0,3) = 1./4.;
  restriction[1](1,3) = 1./4.;
  restriction[1](2,3) = 1./4.;
  restriction[1](3,3) = 1./4.;

  restriction[2](2,2) = 1.0;
  restriction[2](2,1) = 1./2.;
  restriction[2](1,1) = 1./2.;
  restriction[2](2,3) = 1./2.;
  restriction[2](3,3) = 1./2.;
  restriction[2](0,0) = 1./4.;
  restriction[2](1,0) = 1./4.;
  restriction[2](2,0) = 1./4.;
  restriction[2](3,0) = 1./4.;

  restriction[3](3,3) = 1.0;
  restriction[3](0,0) = 1./2.;
  restriction[3](3,0) = 1./2.;
  restriction[3](2,2) = 1./2.;
  restriction[3](3,2) = 1./2.;
  restriction[3](0,1) = 1./4.;
  restriction[3](1,1) = 1./4.;
  restriction[3](2,1) = 1./4.;
  restriction[3](3,1) = 1./4.;

  prolongation[0](0,0) = 1.0;
  prolongation[0](1,0) = 1./2.;
  prolongation[0](1,1) = 1./2.;
  prolongation[0](3,0) = 1./2.;
  prolongation[0](3,3) = 1./2.;
  prolongation[0](2,0) = 1./4.;
  prolongation[0](2,1) = 1./4.;
  prolongation[0](2,2) = 1./4.;
  prolongation[0](2,3) = 1./4.;

  prolongation[1](1,1) = 1.0;
  prolongation[1](0,0) = 1./2.;
  prolongation[1](0,1) = 1./2.;
  prolongation[1](2,1) = 1./2.;
  prolongation[1](2,2) = 1./2.;
  prolongation[1](3,0) = 1./4.;
  prolongation[1](3,1) = 1./4.;
  prolongation[1](3,2) = 1./4.;
  prolongation[1](3,3) = 1./4.;

  prolongation[2](2,2) = 1.0;
  prolongation[2](1,2) = 1./2.;
  prolongation[2](1,1) = 1./2.;
  prolongation[2](3,2) = 1./2.;
  prolongation[2](3,3) = 1./2.;
  prolongation[2](0,0) = 1./4.;
  prolongation[2](0,1) = 1./4.;
  prolongation[2](0,2) = 1./4.;
  prolongation[2](0,3) = 1./4.;

  prolongation[3](3,3) = 1.0;
  prolongation[3](0,0) = 1./2.;
  prolongation[3](0,3) = 1./2.;
  prolongation[3](2,2) = 1./2.;
  prolongation[3](2,3) = 1./2.;
  prolongation[3](1,0) = 1./4.;
  prolongation[3](1,1) = 1./4.;
  prolongation[3](1,2) = 1./4.;
  prolongation[3](1,3) = 1./4.;
};



template <>
inline
double
FELinear<2>::shape_value (const unsigned int i,
			  const Point<2>& p) const
{
  Assert((i<total_dofs), ExcInvalidIndex(i));
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
inline
Point<2>
FELinear<2>::shape_grad (const unsigned int i,
			 const Point<2>& p) const
{
  Assert((i<total_dofs), ExcInvalidIndex(i));
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
void FELinear<2>::get_local_mass_matrix (const DoFHandler<2>::cell_iterator &cell,
					 const Boundary<2> &,
					 dFMatrix &local_mass_matrix) const {
  Assert (local_mass_matrix.n() == total_dofs,
	  ExcWrongFieldDimension(local_mass_matrix.n(),total_dofs));
  Assert (local_mass_matrix.m() == total_dofs,
	  ExcWrongFieldDimension(local_mass_matrix.m(),total_dofs));

/* Get the computation of the local mass matrix by these lines in maple:

   x_real := sum(x[i]*phi[i], i=0..3);
   y_real := sum(y[i]*phi[i], i=0..3);
   phi[0] := (1-xi)*(1-eta);
   phi[1] := xi*(1-eta);
   phi[2] := xi*eta;
   phi[3] := (1-xi)*eta;
   detJ := diff(x_real,xi)*diff(y_real,eta) - diff(x_real,eta)*diff(y_real,xi);

   m := proc (i,j)  int( int(phi[i]*phi[j]*detJ, xi=0..1), eta=0..1); end;

   M := array(0..3,0..3);
   for i from 0 to 3 do
     for j from 0 to 3 do
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
  
  const double t1 = x[1]*y[3],
	       t2 = x[1]*y[2],
	       t3 = x[1]*y[0],
	       t4 = x[0]*y[3],
	       t5 = x[0]*y[1],
	       t6 = x[2]*y[3],
	       t7 = x[3]*y[0],
	       t8 = x[2]*y[1],
	       t9 = x[3]*y[2],
	      t10 = x[3]*y[1],
	      t12 = x[0]*y[2],
	      t13 = x[2]*y[0],
	      t14 = t1/72+t2/36-t3/24-t4/36-t12/72+t5/24+t6/72
		    +t7/36-t8/36-t9/72-t10/72+t13/72,
	      t15 = t2/72-t3/72-t4/72+t5/72+t6/72+t7/72-t8/72-t9/72,
	      t16 = t1/72+t2/72-t3/36-t4/24+t12/72+t5/36+t6/36
		    +t7/24-t8/72-t9/36-t10/72-t13/72,
	      t18 = -t1/72+t2/24-t3/36-t4/72-t12/72+t5/36+t6/36
		    +t7/72-t8/24-t9/36+t10/72+t13/72,
	      t20 = -t1/72+t12/72+t2/36+t5/72-t3/72+t6/24
		    -t9/24-t13/72+t10/72-t4/36+t7/36-t8/36;
  local_mass_matrix(0,0) = t1/18+t2/36-t3/12-t4/12+t5/12+t6/36+t7/12-t8/36-t9/36-t10/18;
  local_mass_matrix(0,1) = t14;
  local_mass_matrix(0,2) = t15;
  local_mass_matrix(0,3) = t16;
  local_mass_matrix(1,0) = t14;
  local_mass_matrix(1,1) = t2/12-t3/12-t4/36-t12/18+t5/12+t6/36+t7/36-t8/12-t9/36+t13/18;
  local_mass_matrix(1,2) = t18;
  local_mass_matrix(1,3) = t15;
  local_mass_matrix(2,0) = t15;
  local_mass_matrix(2,1) = t18;
  local_mass_matrix(2,2) = -t1/18+t2/12+t5/36-t3/36+t6/12-t9/12+t10/18-t4/36+t7/36-t8/12;
  local_mass_matrix(2,3) = t20;
  local_mass_matrix(3,0) = t16;
  local_mass_matrix(3,1) = t15;
  local_mass_matrix(3,2) = t20;
  local_mass_matrix(3,3) = t12/18+t2/36+t5/36-t3/36+t6/12-t9/12-t13/18-t4/12+t7/12-t8/36;
};



template <>
void FELinear<2>::get_face_jacobians (const DoFHandler<2>::face_iterator &face,
				      const Boundary<2>         &,
				      const vector<Point<1> > &unit_points,
				      vector<double>      &face_jacobians) const {
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
void FELinear<2>::get_subface_jacobians (const DoFHandler<2>::face_iterator &face,
					 const unsigned int,
					 const vector<Point<1> > &unit_points,
					 vector<double>      &face_jacobians) const {
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
void FELinear<2>::get_normal_vectors (const DoFHandler<2>::cell_iterator &cell,
				      const unsigned int       face_no,
				      const Boundary<2> &,
				      const vector<Point<1> >  &unit_points,
				      vector<Point<2> >        &normal_vectors) const {
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
void FELinear<2>::get_normal_vectors (const DoFHandler<2>::cell_iterator &cell,
				      const unsigned int       face_no,
				      const unsigned int,
				      const vector<Point<1> >  &unit_points,
				      vector<Point<2> >        &normal_vectors) const {
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
void FELinear<dim>::fill_fe_values (const DoFHandler<dim>::cell_iterator &cell,
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
      for (unsigned int j=0; j<GeometryInfo<dim>::vertices_per_cell; ++j) 
	for (unsigned int l=0; l<n_points; ++l) 
	  q_points[l] += vertices[j] *
			 FELinear<dim>::shape_value(j, unit_points[l]);
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
					       // we want a linear transform and
					       // if we prepend the class name in
					       // front of the #shape_grad#, we
					       // need not use virtual function
					       // calls.
	      const Point<dim> gradient
		= FELinear<dim>::shape_grad (s, unit_points[l]);
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



template <int dim>
void FELinear<dim>::get_ansatz_points (const typename DoFHandler<dim>::cell_iterator &cell,
				       const Boundary<dim>  &,
				       vector<Point<dim> >  &ansatz_points) const {
  Assert (ansatz_points.size() == total_dofs,
	  ExcWrongFieldDimension (ansatz_points.size(), 1<<(dim-1)));
  
  for (unsigned int vertex=0; vertex<GeometryInfo<dim>::vertices_per_cell; ++vertex)
    ansatz_points[vertex] = cell->vertex(vertex);
};



template <int dim>
void FELinear<dim>::get_face_ansatz_points (const typename DoFHandler<dim>::face_iterator &face,
					    const Boundary<dim>  &,
					    vector<Point<dim> >  &ansatz_points) const {
  Assert ((ansatz_points.size() == dofs_per_face) &&
	  (ansatz_points.size() == GeometryInfo<dim>::vertices_per_face),
	  ExcWrongFieldDimension (ansatz_points.size(),
				  GeometryInfo<dim>::vertices_per_face));

  for (unsigned int vertex=0; vertex<dofs_per_face; ++vertex)
    ansatz_points[vertex] = face->vertex(vertex);
};




// explicit instantiations

template class FELinear<deal_II_dimension>;
