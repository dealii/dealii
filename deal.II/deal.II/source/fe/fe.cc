/* $Id$ */

#include <fe/fe.h>
#include <fe/quadrature.h>
#include <grid/tria_iterator.h>
#include <grid/dof_accessor.h>
#include <grid/tria_boundary.h>





/*------------------------------- FiniteElementData ----------------------*/

#if deal_II_dimension == 1

bool FiniteElementData<1>::operator== (const FiniteElementData<1> &f) const {
  return ((dofs_per_vertex == f.dofs_per_vertex) &&
	  (dofs_per_line == f.dofs_per_line) &&
	  (total_dofs == f.total_dofs));
};

#endif



#if deal_II_dimension == 2

bool FiniteElementData<2>::operator== (const FiniteElementData<2> &f) const {
  return ((dofs_per_vertex == f.dofs_per_vertex) &&
	  (dofs_per_line == f.dofs_per_line) &&
	  (dofs_per_quad == f.dofs_per_quad) &&
	  (total_dofs == f.total_dofs));
};

#endif



/*------------------------------- FiniteElementBase ----------------------*/


#if deal_II_dimension == 1

template <>
FiniteElementBase<1>::FiniteElementBase (const unsigned int dofs_per_vertex,
					 const unsigned int dofs_per_line,
					 const unsigned int dofs_per_quad) :
		FiniteElementData<1> (dofs_per_vertex,
				      dofs_per_line)
{
  Assert (dofs_per_quad==0, ExcInternalError());

  const unsigned int dim=1;
  for (unsigned int i=0; i<GeometryInfo<dim>::children_per_cell; ++i) 
    {
      restriction[i].reinit (total_dofs, total_dofs);
      prolongation[i].reinit (total_dofs, total_dofs);
    };
  interface_constraints.reinit (1,1);
  interface_constraints(0,0)=1.;
};

#endif


#if deal_II_dimension == 2

template <>
FiniteElementBase<2>::FiniteElementBase (const unsigned int dofs_per_vertex,
					 const unsigned int dofs_per_line,
					 const unsigned int dofs_per_quad) :
		FiniteElementData<2> (dofs_per_vertex,
				      dofs_per_line,
				      dofs_per_quad)
{
  const unsigned int dim=2;
  for (unsigned int i=0; i<GeometryInfo<dim>::children_per_cell; ++i) 
    {
      restriction[i].reinit (total_dofs, total_dofs);
      prolongation[i].reinit (total_dofs, total_dofs);
    };
  interface_constraints.reinit (dofs_per_vertex+2*dofs_per_line,
				2*dofs_per_vertex+dofs_per_line);
};

#endif



template <int dim>
const dFMatrix &
FiniteElementBase<dim>::restrict (const unsigned int child) const {
  Assert (child<GeometryInfo<dim>::children_per_cell, ExcInvalidIndex(child));
  return restriction[child];
};



template <int dim>
const dFMatrix &
FiniteElementBase<dim>::prolongate (const unsigned int child) const {
  Assert (child<GeometryInfo<dim>::children_per_cell, ExcInvalidIndex(child));
  return prolongation[child];
};



template <int dim>
const dFMatrix &
FiniteElementBase<dim>::constraints () const {
  if (dim==1)
    Assert ((interface_constraints.m()==1) && (interface_constraints.n()==1),
	    ExcWrongInterfaceMatrixSize(interface_constraints.m(),
					interface_constraints.n()));
  
  return interface_constraints;
};



template <int dim>
bool FiniteElementBase<dim>::operator == (const FiniteElementBase<dim> &f) const {
  return ((static_cast<FiniteElementData<dim> >(*this) ==
	   static_cast<FiniteElementData<dim> >(f)) &&
	  (interface_constraints == f.interface_constraints));
};





/*------------------------------- FiniteElement ----------------------*/

// declare this function to be explicitely specialized before first use
// egcs wants this, but gcc2.8.1 produces an internal compiler error, so
// we drop this declaration again for the time being

#if deal_II_dimension == 1

//template <>
//void FiniteElement<1>::get_ansatz_points (const DoFHandler<1>::cell_iterator &cell,
//					  const Boundary<1> &,
//					  vector<Point<1> > &ansatz_points) const;


template <>
void FiniteElement<1>::fill_fe_values (const DoFHandler<1>::cell_iterator &cell,
				       const vector<Point<1> > &unit_points,
				       vector<dFMatrix>  &jacobians,
				       const bool         compute_jacobians,
				       vector<Point<1> > &ansatz_points,
				       const bool         compute_ansatz_points,
				       vector<Point<1> > &q_points,
				       const bool         compute_q_points,
				       const Boundary<1> &boundary) const {
  Assert (jacobians.size() == unit_points.size(),
	  ExcWrongFieldDimension(jacobians.size(), unit_points.size()));
  Assert (q_points.size() == unit_points.size(),
	  ExcWrongFieldDimension(q_points.size(), unit_points.size()));
  Assert (ansatz_points.size() == total_dofs,
	  ExcWrongFieldDimension(ansatz_points.size(), total_dofs));


				   // local mesh width
  const double h=(cell->vertex(1)(0) - cell->vertex(0)(0));

  for (unsigned int i=0; i<q_points.size(); ++i) 
    {
      if (compute_jacobians)
	jacobians[i](0,0) = 1./h;
      if (compute_q_points)
	q_points[i] = cell->vertex(0) + h*unit_points[i];
    };

				   // compute ansatz points. The first ones
				   // belong to vertex one, the second ones
				   // to vertex two, all following are
				   // equally spaced along the line
  if (compute_ansatz_points)
    get_ansatz_points (cell, boundary, ansatz_points);
};



template <>
void FiniteElement<1>::fill_fe_face_values (const DoFHandler<1>::cell_iterator &,
					    const unsigned int       ,
					    const vector<Point<0> > &,
					    const vector<Point<1> > &,
					    vector<dFMatrix>        &,
					    const bool               ,
					    vector<Point<1> >       &,
					    const bool               ,
					    vector<Point<1> >       &,
					    const bool               ,
					    vector<double>          &,
					    const bool              ,
					    vector<Point<1> >       &,
					    const bool,
					    const Boundary<1>       &) const {
  Assert (false, ExcNotImplemented());
};


template <>
void FiniteElement<1>::fill_fe_subface_values (const DoFHandler<1>::cell_iterator &,
					       const unsigned int       ,
					       const unsigned int       ,
					       const vector<Point<0> > &,
					       const vector<Point<1> > &,
					       vector<dFMatrix>        &,
					       const bool               ,
					       vector<Point<1> >       &,
					       const bool               ,
					       vector<double>          &,
					       const bool               ,
					       vector<Point<1> >       &,
					       const bool,
					       const Boundary<1>       &) const {
  Assert (false, ExcNotImplemented());
};



template <>
void FiniteElement<1>::get_ansatz_points (const DoFHandler<1>::cell_iterator &cell,
					  const Boundary<1> &,
					  vector<Point<1> > &ansatz_points) const {
  Assert (ansatz_points.size() == total_dofs,
	  ExcWrongFieldDimension(ansatz_points.size(), total_dofs));
				   // compute ansatz points. The first ones
				   // belong to vertex one, the second ones
				   // to vertex two, all following are
				   // equally spaced along the line
  unsigned int next = 0;
				   // local mesh width
  const double h=(cell->vertex(1)(0) - cell->vertex(0)(0));
				   // first the dofs in the vertices
  for (unsigned int vertex=0; vertex<2; vertex++) 
    for (unsigned int i=0; i<dofs_per_vertex; ++i)
      ansatz_points[next++] = cell->vertex(vertex);
  
				   // now dofs on line
  for (unsigned int i=0; i<dofs_per_line; ++i) 
    ansatz_points[next++] = cell->vertex(0) +
			    Point<1>((i+1.0)/(total_dofs+1.0)*h);
};

#endif



template <int dim>
void FiniteElement<dim>::fill_fe_values (const DoFHandler<dim>::cell_iterator &,
					 const vector<Point<dim> > &,
					 vector<dFMatrix> &,
					 const bool,
					 vector<Point<dim> > &,
					 const bool,
					 vector<Point<dim> > &,
					 const bool,
					 const Boundary<dim> &) const {
  Assert (false, ExcPureFunctionCalled());
};



template <int dim>
void FiniteElement<dim>::fill_fe_face_values (const DoFHandler<dim>::cell_iterator &cell,
					      const unsigned int           face_no,
					      const vector<Point<dim-1> > &unit_points,
					      const vector<Point<dim> > &global_unit_points,
					      vector<dFMatrix>    &jacobians,
					      const bool           compute_jacobians,
					      vector<Point<dim> > &ansatz_points,
					      const bool           compute_ansatz_points,
					      vector<Point<dim> > &q_points,
					      const bool           compute_q_points,
					      vector<double>      &face_jacobi_determinants,
					      const bool           compute_face_jacobians,
					      vector<Point<dim> > &normal_vectors,
					      const bool           compute_normal_vectors,
					      const Boundary<dim> &boundary) const {
  Assert (jacobians.size() == unit_points.size(),
	  ExcWrongFieldDimension(jacobians.size(), unit_points.size()));
  Assert (q_points.size() == unit_points.size(),
	  ExcWrongFieldDimension(q_points.size(), unit_points.size()));
  Assert (global_unit_points.size() == unit_points.size(),
	  ExcWrongFieldDimension(global_unit_points.size(), unit_points.size()));
  Assert (ansatz_points.size() == dofs_per_face,
	  ExcWrongFieldDimension(ansatz_points.size(), dofs_per_face));
  
  static vector<Point<dim> > dummy(total_dofs);
  fill_fe_values (cell, global_unit_points,
		  jacobians, compute_jacobians,
		  dummy, false,
		  q_points, compute_q_points,
		  boundary);
  
  if (compute_ansatz_points)
    get_face_ansatz_points (cell->face(face_no), boundary, ansatz_points);

  if (compute_face_jacobians)
    get_face_jacobians (cell->face(face_no), boundary,
			unit_points, face_jacobi_determinants);

  if (compute_normal_vectors)
    get_normal_vectors (cell, face_no, boundary,
			unit_points, normal_vectors);
};




template <int dim>
void FiniteElement<dim>::fill_fe_subface_values (const DoFHandler<dim>::cell_iterator &cell,
						 const unsigned int           face_no,
						 const unsigned int           subface_no,
						 const vector<Point<dim-1> > &unit_points,
						 const vector<Point<dim> > &global_unit_points,
						 vector<dFMatrix>    &jacobians,
						 const bool           compute_jacobians,
						 vector<Point<dim> > &q_points,
						 const bool           compute_q_points,
						 vector<double>      &face_jacobi_determinants,
						 const bool           compute_face_jacobians,
						 vector<Point<dim> > &normal_vectors,
						 const bool           compute_normal_vectors,
						 const Boundary<dim> &boundary) const {
  Assert (jacobians.size() == unit_points.size(),
	  ExcWrongFieldDimension(jacobians.size(), unit_points.size()));
  Assert (q_points.size() == unit_points.size(),
	  ExcWrongFieldDimension(q_points.size(), unit_points.size()));
  Assert (global_unit_points.size() == unit_points.size(),
	  ExcWrongFieldDimension(global_unit_points.size(), unit_points.size()));

  static vector<Point<dim> > dummy(total_dofs);
  fill_fe_values (cell, global_unit_points,
		  jacobians, compute_jacobians,
		  dummy, false,
		  q_points, compute_q_points,
		  boundary);
  
  if (compute_face_jacobians)
    get_subface_jacobians (cell->face(face_no), subface_no,
			   unit_points, face_jacobi_determinants);

  if (compute_normal_vectors)
    get_normal_vectors (cell, face_no, subface_no,
			unit_points, normal_vectors);
};



template <int dim>
void FiniteElement<dim>::get_ansatz_points (const DoFHandler<dim>::cell_iterator &,
					    const Boundary<dim> &,
					    vector<Point<dim> > &) const {
  Assert (false, ExcPureFunctionCalled());
};



#if deal_II_dimension == 1

template <>
inline
double
FELinearMapping<1>::linear_shape_value(const unsigned int i,
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
Point<1>
FELinearMapping<1>::linear_shape_grad(const unsigned int i,
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
void FEQuadraticSub<1>::get_face_jacobians (const DoFHandler<1>::face_iterator &,
					    const Boundary<1>         &,
					    const vector<Point<0> > &,
					    vector<double>      &) const {
  Assert (false, ExcInternalError());
};



template <>
void FEQuadraticSub<1>::get_subface_jacobians (const DoFHandler<1>::face_iterator &,
					       const unsigned int           ,
					       const vector<Point<0> > &,
					       vector<double>      &) const {
  Assert (false, ExcInternalError());
};



template <>
void FEQuadraticSub<1>::get_normal_vectors (const DoFHandler<1>::cell_iterator &,
					    const unsigned int,
					    const Boundary<1> &,
					    const vector<Point<0> > &,
					    vector<Point<1> > &) const {
  Assert (false, ExcInternalError());
};



template <>
void FEQuadraticSub<1>::get_normal_vectors (const DoFHandler<1>::cell_iterator &,
					    const unsigned int,
					    const unsigned int,
					    const vector<Point<0> > &,
					    vector<Point<1> > &) const {
  Assert (false, ExcInternalError());
};


template <>
void FELinearMapping<1>::fill_fe_values (const DoFHandler<1>::cell_iterator &cell,
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



#endif



#if deal_II_dimension == 2

template <>
inline
double
FELinearMapping<2>::linear_shape_value (const unsigned int i,
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
Point<2>
FELinearMapping<2>::linear_shape_grad (const unsigned int i,
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



template <int dim>
void FELinearMapping<dim>::fill_fe_values (const DoFHandler<dim>::cell_iterator &cell,
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



/*------------------------------- Explicit Instantiations -------------*/

template class FiniteElementData<deal_II_dimension>;
template class FiniteElementBase<deal_II_dimension>;
template class FiniteElement<deal_II_dimension>;
template class FELinearMapping<deal_II_dimension>;


