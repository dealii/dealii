/* $Id$ */
/* Copyright W. Bangerth, University of Heidelberg, 1998 */


#include <fe/fe_linear_mapping.h>
#include <base/quadrature.h>
#include <grid/tria_iterator.h>
#include <grid/dof_accessor.h>
#include <grid/tria_boundary.h>

#include <cmath>



/*---------------------------- FELinearMapping ----------------------------------*/



#if deal_II_dimension == 1

template <>
FELinearMapping<1>::FELinearMapping (const unsigned int dofs_per_vertex,
				     const unsigned int dofs_per_line,
				     const unsigned int dofs_per_quad,
				     const unsigned int dofs_per_hex,
				     const unsigned int n_components) :
		FiniteElement<1> (FiniteElementData<1> (dofs_per_vertex,
							dofs_per_line,
							GeometryInfo<1>::vertices_per_cell,
							n_components))
{
  Assert (dofs_per_quad==0, ExcInvalidData());
  Assert (dofs_per_hex==0,  ExcInvalidData());
};



template <>
inline
double
FELinearMapping<1>::shape_value_transform (const unsigned int i,
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
inline
Tensor<1,1>
FELinearMapping<1>::shape_grad_transform(const unsigned int i,
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
void FELinearMapping<1>::get_face_jacobians (const DoFHandler<1>::face_iterator &,
					     const Boundary<1>         &,
					     const vector<Point<0> > &,
					     vector<double>      &) const {
  Assert (false, ExcInternalError());
};



template <>
void FELinearMapping<1>::get_subface_jacobians (const DoFHandler<1>::face_iterator &,
						const unsigned int           ,
						const vector<Point<0> > &,
						vector<double>      &) const {
  Assert (false, ExcInternalError());
};



template <>
void FELinearMapping<1>::get_normal_vectors (const DoFHandler<1>::cell_iterator &,
					     const unsigned int,
					     const Boundary<1> &,
					     const vector<Point<0> > &,
					     vector<Point<1> > &) const {
  Assert (false, ExcInternalError());
};



template <>
void FELinearMapping<1>::get_normal_vectors (const DoFHandler<1>::cell_iterator &,
					    const unsigned int,
					    const unsigned int,
					    const vector<Point<0> > &,
					    vector<Point<1> > &) const {
  Assert (false, ExcInternalError());
};


template <>
void FELinearMapping<1>::fill_fe_values (const DoFHandler<1>::cell_iterator &cell,
					 const vector<Point<1> > &unit_points,
					 vector<Tensor<2,1> >    &jacobians,
					 const bool            compute_jacobians,
					 vector<Tensor<3,1> > &jacobians_grad,
					 const bool            compute_jacobians_grad,
					 vector<Point<1> >    &support_points,
					 const bool            compute_support_points,
					 vector<Point<1> >    &q_points,
					 const bool            compute_q_points,
					 const FullMatrix<double>       &shape_values_transform,
					 const vector<vector<Tensor<1,1> > > &shape_gradients_transform,
					 const Boundary<1> &boundary) const {
				   // simply pass down
  FiniteElement<1>::fill_fe_values (cell, unit_points,
				    jacobians, compute_jacobians,
				    jacobians_grad, compute_jacobians_grad,
				    support_points, compute_support_points,
				    q_points, compute_q_points,
				    shape_values_transform, shape_gradients_transform,
				    boundary);
};



#endif



#if deal_II_dimension == 2

template <>
FELinearMapping<2>::FELinearMapping (const unsigned int dofs_per_vertex,
				     const unsigned int dofs_per_line,
				     const unsigned int dofs_per_quad,
				     const unsigned int dofs_per_hex,
				     const unsigned int n_components) :
		FiniteElement<2> (FiniteElementData<2> (dofs_per_vertex,
							dofs_per_line,
							dofs_per_quad,
							GeometryInfo<2>::vertices_per_cell,
							n_components))
{
  Assert (dofs_per_hex == 0, ExcInvalidData());
};



template <>
inline
double
FELinearMapping<2>::shape_value_transform (const unsigned int i,
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
inline
Tensor<1,2>
FELinearMapping<2>::shape_grad_transform (const unsigned int i,
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
void FELinearMapping<2>::get_face_jacobians (const DoFHandler<2>::face_iterator &face,
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
void FELinearMapping<2>::get_subface_jacobians (const DoFHandler<2>::face_iterator &face,
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
void FELinearMapping<2>::get_normal_vectors (const DoFHandler<2>::cell_iterator &cell,
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
void FELinearMapping<2>::get_normal_vectors (const DoFHandler<2>::cell_iterator &cell,
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



#if deal_II_dimension == 3

template <>
FELinearMapping<3>::FELinearMapping (const unsigned int dofs_per_vertex,
				     const unsigned int dofs_per_line,
				     const unsigned int dofs_per_quad,
				     const unsigned int dofs_per_hex,
				     const unsigned int n_components) :
		FiniteElement<3> (FiniteElementData<3> (dofs_per_vertex,
							dofs_per_line,
							dofs_per_quad,
							dofs_per_hex,
							GeometryInfo<3>::vertices_per_cell,
							n_components))
{};



template <>
inline
double
FELinearMapping<3>::shape_value_transform (const unsigned int i,
					   const Point<3>& p) const
{
  Assert((i<8), ExcInvalidIndex(i));
  switch (i)
    {
      case 0: return 1.0-p(0)+(-1.0+p(0))*p(1)+(-1.0+p(0)+(1.0-p(0))*p(1))*p(2);
      case 1: return p(0)-p(0)*p(1)+(-p(0)+p(0)*p(1))*p(2);
      case 2: return (p(0)-p(0)*p(1))*p(2);
      case 3: return (1.0-p(0)+(-1.0+p(0))*p(1))*p(2);
      case 4: return (1.0-p(0))*p(1)+(-1.0+p(0))*p(1)*p(2);
      case 5: return p(0)*p(1)-p(0)*p(1)*p(2);
      case 6: return p(0)*p(1)*p(2);
      case 7: return (1.0-p(0))*p(1)*p(2);
    }
  return 0.;
};



template <>
inline
Tensor<1,3>
FELinearMapping<3>::shape_grad_transform (const unsigned int i,
					  const Point<3>& p) const
{
  Assert((i<8), ExcInvalidIndex(i));
  switch (i)
    {
      case 0: return Point<3>(-1.0+p(1)+(1.0-p(1))*p(2),
			      -1.0+p(0)+(1.0-p(0))*p(2),
			      -1.0+p(0)+(1.0-p(0))*p(1));
      case 1: return Point<3>(1.0-p(1)+(-1.0+p(1))*p(2),
			      -p(0)+p(0)*p(2),
			      -p(0)+p(0)*p(1));
      case 2: return Point<3>((1.0-p(1))*p(2),
			      -p(0)*p(2),
			      p(0)-p(0)*p(1));
      case 3: return Point<3>((-1.0+p(1))*p(2),
			      (-1.0+p(0))*p(2),
			      1.0-p(0)+(-1.0+p(0))*p(1));
      case 4: return Point<3>(-p(1)+p(1)*p(2),
			      1.0-p(0)+(-1.0+p(0))*p(2),
			      (-1.0+p(0))*p(1));
      case 5: return Point<3>(p(1)-p(1)*p(2),
			      p(0)-p(0)*p(2),
			      -p(0)*p(1));
      case 6: return Point<3>(p(1)*p(2),
			      p(0)*p(2),
			      p(0)*p(1));
      case 7: return Point<3>(-p(1)*p(2),
			      (1.0-p(0))*p(2),
			      (1.0-p(0))*p(1));
    }
  return Point<3> ();
};



template <>
void FELinearMapping<3>::get_face_jacobians (const DoFHandler<3>::face_iterator &face,
					     const Boundary<3>         &,
					     const vector<Point<2> > &unit_points,
					     vector<double> &face_jacobians) const {
  Assert (unit_points.size() == face_jacobians.size(),
	  ExcWrongFieldDimension (unit_points.size(), face_jacobians.size()));

				   // the computation of the face jacobians is
				   // along the following lines: let x_i be
				   // the four vertices of a face, then the
				   // unit point (xi,eta) is mapped to the
				   // point vec x(xi,eta)=\sum x_i phi_i(xi,eta),
				   // with phi_i being the shape functions
				   // of this face
				   //
				   // now, while d(xi) d(eta) is the area
				   // element on the unit face,
				   // abs(dx dy) is the respective element
				   // of the real face. to compute it, we
				   // compute the image of the elements d(xi)
				   // and d(eta):
				   // (\vec x(xi+dxi,eta) - \vec x(xi,eta) )
  				   // (\vec x(xi,eta+deta) - \vec x(xi,eta) )
				   // the area then is the norm of the
				   // cross product of these two vectors
				   // and the determinant is the area
				   // divided by d(xi)d(eta)
				   //
				   // written down, we remark that the
				   // determinant we are looking for is
				   // the cross product of the following
				   // two vectors:
				   // d/d(xi)  vec x(xi,eta)
				   // d/d(eta) vec x(xi,eta)
				   // we then arrive at:
				   //
				   // detJ =
				   // || \sum_l \sum_k \phi_{l,xi} \phi_{k,eta}
				   //        x_l \times x_k ||
				   //
				   // a maple script doing this computation is
				   // in the <scripts> directory
  const Point<3> vertices[4] = { face->vertex(0),
				 face->vertex(1),
				 face->vertex(2),
				 face->vertex(3)   };

  for (unsigned int point=0; point<unit_points.size(); ++point)
    {
      const double xi  = unit_points[point](0),
		   eta = unit_points[point](1);

      const double t1 = 1.0-eta;
      const double t6 = -vertices[0](1)*t1+vertices[1](1)*t1+vertices[2](1)*eta-vertices[3](1)*eta;
      const double t7 = 1.0-xi;
      const double t12 = -vertices[0](2)*t7-vertices[1](2)*xi+vertices[2](2)*xi+vertices[3](2)*t7;
      const double t18 = -vertices[0](2)*t1+vertices[1](2)*t1+vertices[2](2)*eta-vertices[3](2)*eta;
      const double t23 = -vertices[0](1)*t7-vertices[1](1)*xi+vertices[2](1)*xi+vertices[3](1)*t7;
      const double t26 = fabs(t6*t12-t18*t23);
      const double t27 = t26*t26;
      const double t32 = -vertices[0](0)*t7-vertices[1](0)*xi+vertices[2](0)*xi+vertices[3](0)*t7;
      const double t38 = -vertices[0](0)*t1+vertices[1](0)*t1+vertices[2](0)*eta-vertices[3](0)*eta;
      const double t41 = fabs(t18*t32-t38*t12);
      const double t42 = t41*t41;
      const double t46 = fabs(t38*t23-t6*t32);
      const double t47 = t46*t46;
      face_jacobians[point] = sqrt(t27+t42+t47);
    };
};



template <>
void FELinearMapping<3>::get_subface_jacobians (const DoFHandler<3>::face_iterator &/*face*/,
						const unsigned int           ,
						const vector<Point<2> > &unit_points,
						vector<double> &face_jacobians) const {
  Assert (false,
	  ExcWrongFieldDimension (unit_points.size(), face_jacobians.size()));
};



template <>
void FELinearMapping<3>::get_normal_vectors (const DoFHandler<3>::cell_iterator &cell,
					     const unsigned int       face_no,
					     const Boundary<3>       &,
					     const vector<Point<2> > &unit_points,
					     vector<Point<3> > &normal_vectors) const {
  Assert (unit_points.size() == normal_vectors.size(),
	  ExcWrongFieldDimension (unit_points.size(), normal_vectors.size()));
  
				   // taken from the same script as is
				   // the computation of the jacobian
				   // determinant above

  const Point<3> vertices[4] = { cell->face(face_no)->vertex(0),
				 cell->face(face_no)->vertex(1),
				 cell->face(face_no)->vertex(2),
				 cell->face(face_no)->vertex(3)   };

  for (unsigned int point=0; point<unit_points.size(); ++point)
    {
      const double xi  = unit_points[point](0),
		   eta = unit_points[point](1);
      
      const double t1 = 1.0-eta;
      const double t6 = -vertices[0](1)*t1+vertices[1](1)*t1+vertices[2](1)*eta-vertices[3](1)*eta;
      const double t7 = 1.0-xi;
      const double t12 = -vertices[0](2)*t7-vertices[1](2)*xi+vertices[2](2)*xi+vertices[3](2)*t7;
      const double t18 = -vertices[0](2)*t1+vertices[1](2)*t1+vertices[2](2)*eta-vertices[3](2)*eta;
      const double t23 = -vertices[0](1)*t7-vertices[1](1)*xi+vertices[2](1)*xi+vertices[3](1)*t7;
      const double t25 = t6*t12-t18*t23;
      const double t26 = fabs(t25);
      const double t27 = t26*t26;
      const double t32 = -vertices[0](0)*t7-vertices[1](0)*xi+vertices[2](0)*xi+vertices[3](0)*t7;
      const double t38 = -vertices[0](0)*t1+vertices[1](0)*t1+vertices[2](0)*eta-vertices[3](0)*eta;
      const double t40 = t18*t32-t38*t12;
      const double t41 = fabs(t40);
      const double t42 = t41*t41;
      const double t45 = t38*t23-t6*t32;
      const double t46 = fabs(t45);
      const double t47 = t46*t46;
      const double t49 = sqrt(t27+t42+t47);
      const double t50 = 1/t49;
      
      normal_vectors[point](0) = t25*t50;
      normal_vectors[point](1) = t40*t50;
      normal_vectors[point](2) = t45*t50;

      if ((face_no == 1) ||
	  (face_no == 2) ||
	  (face_no == 5))
	normal_vectors[point] *= -1;
    };  
};



template <>
void FELinearMapping<3>::get_normal_vectors (const DoFHandler<3>::cell_iterator &/*cell*/,
					     const unsigned int       /*face_no*/,
					     const unsigned int,
					     const vector<Point<2> > &unit_points,
					     vector<Point<3> > &normal_vectors) const {
				   // more or less copied from the linear
				   // finite element
				   // note, that in 2D the normal vectors to the
				   // subface have the same direction as that
				   // for the face
  Assert (false,
	  ExcWrongFieldDimension (unit_points.size(), normal_vectors.size()));
};

#endif



template <int dim>
void FELinearMapping<dim>::fill_fe_values (const DoFHandler<dim>::cell_iterator &cell,
					   const vector<Point<dim> >            &unit_points,
					   vector<Tensor<2,dim> >               &jacobians,
					   const bool              compute_jacobians,
					   vector<Tensor<3,dim> > &jacobians_grad,
					   const bool              compute_jacobians_grad,
					   vector<Point<dim> > &support_points,
					   const bool           compute_support_points,
					   vector<Point<dim> > &q_points,
					   const bool           compute_q_points,
					   const FullMatrix<double>      &shape_values_transform,
					   const vector<vector<Tensor<1,dim> > > &/*shape_grad_transform*/,
					   const Boundary<dim> &boundary) const
{
  Assert ((!compute_jacobians) || (jacobians.size() == unit_points.size()),
	  ExcWrongFieldDimension(jacobians.size(), unit_points.size()));
  Assert ((!compute_jacobians_grad) || (jacobians_grad.size() == unit_points.size()),
	  ExcWrongFieldDimension(jacobians_grad.size(), unit_points.size()));
  Assert ((!compute_q_points) || (q_points.size() == unit_points.size()),
	  ExcWrongFieldDimension(q_points.size(), unit_points.size()));
  Assert ((!compute_support_points) || (support_points.size() == total_dofs),
	  ExcWrongFieldDimension(support_points.size(), total_dofs));

  
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
	  q_points[l] += vertices[j] * shape_values_transform(j, l);
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

  Indeed, this was the old way we did it; the code is below. However, there
  is a more efficient way, namely setting up M analytically, doing the
  inversion analyically and only then doing the evaluation; a small Maple
  script does this (it is part of the <lagrange> script in the <scripts>
  subdirectory).
  
  if (compute_jacobians) 
    {
      FullMatrix<double> M(dim,dim);
      for (unsigned int l=0; l<n_points; ++l) 
	{
	  M.clear ();
	  for (unsigned int s=0; s<GeometryInfo<dim>::vertices_per_cell; ++s)
	    {
					       // we want the linear transform,
					       // so use that function
	      const Point<dim> gradient = shape_grad_transform[s][l];
	      for (unsigned int i=0; i<dim; ++i)
		for (unsigned int j=0; j<dim; ++j)
		  M(i,j) += vertices[s](i) * gradient(j);
	    };
	  jacobians[l].invert(M);
	};
    };

    
  One last note regarding whether we have to invert M or M transposed: it is
  easy to try out, by computing the gradients of a function on a distorted
  cell (just move one vertex) where the nodal values for linear elements
  are one for the moved vertex and zero otherwise. Please also note that
  when computing the gradients on the real cell, the jacobian matrix
  is multiplied to the unit cell gradient *from the right*! be very careful
  with these things.

  The following little program tests the correct behaviour. The program can
  also be found in the <tests> directory.

  -------------------------------------------
  #include <grid/tria.h>
  #include <grid/tria_boundary.h>
  #include <grid/dof.h>
  #include <fe/fe_values.h>
  #include <fe/fe_lib.lagrange.h>
  #include <base/quadrature_lib.h>
  #include <grid/tria_iterator.h>
  #include <grid/dof_accessor.h>
  #include <lac/vector.h>

  int main () {
    Triangulation<2> tria;
    tria.create_hypercube (0,1);
    tria.begin_active()->vertex(2)(0) = 2;

    DoFHandler<2> dof(&tria);
    FELinear<2> fe;
    dof.distribute_dofs(fe);

    StraightBoundary<2> b;
    QTrapez<2> q;
    FEValues<2> fevalues(fe,q,update_gradients);
    fevalues.reinit (dof.begin_active(),b);
  
  
    Vector<double> val(4);
    val(2) = 1;

    vector<Point<2> > grads(4);
    fevalues.get_function_grads (val, grads);

    for (unsigned int i=0; i<4; ++i)
      cout << "Vertex " << i
	   << "   grad=" << grads[i] << endl;
  };
  ---------------------------------------------
  
  The correct output should be
  --------------------------------
  Vertex 0   grad=0 0
  Vertex 1   grad=0.5 0
  Vertex 2   grad=0 1
  Vertex 3   grad=0.5 0.5
  --------------------------------
  and the wrong would be
  --------------------------------
  Vertex 0   grad=0 0
  Vertex 1   grad=0.5 0
  Vertex 2   grad=-1 1
  Vertex 3   grad=0 1
  --------------------------------  
*/

  if (compute_jacobians)
    compute_jacobian_matrices (cell, unit_points, jacobians);
  
  if (compute_jacobians_grad)
    compute_jacobian_gradients (cell, unit_points, jacobians_grad);
    
  if (compute_support_points)
    get_support_points (cell, boundary, support_points);
};




/*------------------------------- Explicit Instantiations -------------*/

template class FELinearMapping<deal_II_dimension>;
