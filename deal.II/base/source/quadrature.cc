/* $Id$ */
/* Copyright W. Bangerth, University of Heidelberg, 1998 */


#include <base/quadrature.h>


template <int dim>
Quadrature<dim>::Quadrature (const unsigned int n_q) :
		n_quadrature_points(n_q),
		quadrature_points (n_q, Point<dim>()),
		weights (n_q, 0) {};


template <>
Quadrature<1>::Quadrature (const Quadrature<0> &,
			   const Quadrature<1> &) {
  Assert (false, ExcInternalError());
};



template <int dim>
Quadrature<dim>::Quadrature (const Quadrature<dim-1> &q1,
			     const Quadrature<1>     &q2) :
		n_quadrature_points (q1.n_quadrature_points *
				     q2.n_quadrature_points),
		quadrature_points (n_quadrature_points),
		weights (n_quadrature_points, 0)
{
  unsigned int present_index = 0;
  for (unsigned int i=0; i<q1.n_quadrature_points; ++i)
    for (unsigned int j=0; j<q2.n_quadrature_points; ++j)
      {
					 // compose coordinates of
					 // new quadrature point by tensor
					 // product in the last component
	for (unsigned int d=0; d<dim-1; ++d)
	  quadrature_points[present_index](d)
	    = q1.quad_point(i)(d);
	quadrature_points[present_index](dim-1)
	  = q2.quad_point(j)(0);
					       
	weights[present_index] = q1.weight(i) * q2.weight(j);

	++present_index;
      };

#ifdef DEBUG
  double sum = 0;
  for (unsigned int i=0; i<n_quadrature_points; ++i)
    sum += weights[i];
				   // we cant guarantee the sum of weights
				   // to be exactly one, but it should be
				   // near that. 
  Assert ((sum>0.999999) && (sum<1.000001), ExcInternalError());
#endif
};



template <int dim>
Quadrature<dim>::~Quadrature () {};



template <int dim>
const Point<dim> & Quadrature<dim>::quad_point (const unsigned int i) const {
  Assert (i<n_quadrature_points, ExcInvalidIndex(i, n_quadrature_points));
  return quadrature_points[i];
};



template <int dim>
const vector<Point<dim> > & Quadrature<dim>::get_quad_points () const {
  return quadrature_points;
};



template <int dim>
double Quadrature<dim>::weight (const unsigned int i) const {
  Assert (i<n_quadrature_points, ExcInvalidIndex(i, n_quadrature_points));
  return weights[i];
};



template <int dim>
const vector<double> & Quadrature<dim>::get_weights () const {
  return weights;
};



template <>
void QProjector<2>::project_to_face (const Quadrature<1> &quadrature,
				     const unsigned int   face_no,
				     vector<Point<2> >   &q_points) {
  const unsigned int dim=2;
  Assert (face_no<2*dim, ExcInvalidIndex (face_no, 2*dim));
  
  for (unsigned int p=0; p<quadrature.n_quadrature_points; ++p)
    switch (face_no)
      {
	case 0:
	      q_points[p] = Point<dim>(quadrature.quad_point(p)(0),0);
	      break;	   
	case 1:
	      q_points[p] = Point<dim>(1,quadrature.quad_point(p)(0));
	      break;	   
	case 2:
	      q_points[p] = Point<dim>(quadrature.quad_point(p)(0),1);
	      break;	   
	case 3:
	      q_points[p] = Point<dim>(0,quadrature.quad_point(p)(0));
	      break;
	default:
	      Assert (false, ExcInternalError());
      };
};



template <>
void QProjector<2>::project_to_subface (const Quadrature<1> &quadrature,
					const unsigned int   face_no,
					const unsigned int   subface_no,
					vector<Point<2> >   &q_points) {
  const unsigned int dim=2;
  Assert (face_no<2*dim, ExcInvalidIndex (face_no, 2*dim));
  Assert (subface_no<(1<<(dim-1)), ExcInvalidIndex (face_no, 1<<(dim-1)));
  
  for (unsigned int p=0; p<quadrature.n_quadrature_points; ++p)
    switch (face_no)
      {
	case 0:
	      switch (subface_no) 
		{
		  case 0:
			q_points[p]
			  = Point<dim>(quadrature.quad_point(p)(0)/2,0);
			break;
		  case 1:
			q_points[p]
			  = Point<dim>(quadrature.quad_point(p)(0)/2+0.5,0);
			break;
		  default:
			Assert (false, ExcInternalError());
		};
	      break;	   
	case 1:
	      switch (subface_no) 
		{
		  case 0:
			q_points[p] = Point<dim>(1,quadrature.quad_point(p)(0)/2);
			break;
		  case 1:
			q_points[p] = Point<dim>(1,quadrature.quad_point(p)(0)/2+0.5);
			break;
		  default:
			Assert (false, ExcInternalError());
		};
	      break;	   
	case 2:
	      switch (subface_no) 
		{
		  case 0:
			q_points[p] = Point<dim>(quadrature.quad_point(p)(0)/2,1);
			break;
		  case 1:
			q_points[p] = Point<dim>(quadrature.quad_point(p)(0)/2+0.5,1);
			break;
		  default:
			Assert (false, ExcInternalError());
		};
	      break;	   
	case 3:
	      switch (subface_no) 
		{
		  case 0:
			q_points[p] = Point<dim>(0,quadrature.quad_point(p)(0)/2);
			break;
		  case 1:
			q_points[p] = Point<dim>(0,quadrature.quad_point(p)(0)/2+0.5);
			break;
		  default:
			Assert (false, ExcInternalError());
		};
	      break;
	default:
	      Assert (false, ExcInternalError());
      };
};




// explicit instantiations; note: we need them all for all dimensions
template class Quadrature<1>;
template class Quadrature<2>;
