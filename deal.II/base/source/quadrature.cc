//----------------------------  quadrature.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  quadrature.cc  ---------------------------


#include <grid/geometry_info.h>
#include <base/quadrature.h>
#include <base/memory_consumption.h>
#include <cmath>


template <>
Quadrature<0>::Quadrature (const unsigned int)
              : n_quadrature_points(0)
{};



template <int dim>
Quadrature<dim>::Quadrature (const unsigned int n_q) :
		n_quadrature_points(n_q),
		quadrature_points (n_q, Point<dim>()),
		weights (n_q, 0)
{}



// if necessary try to work around a bug in the IBM xlC compiler
#ifdef XLC_WORK_AROUND_STD_BUG
using namespace std;
#endif


template <int dim>
Quadrature<dim>::Quadrature (const typename std::vector<Point<dim> > &points,
			     const std::vector<double>               &weights)
		:
		n_quadrature_points(points.size()),
		quadrature_points(points),
		weights(weights)
{
  Assert(weights.size() == points.size(),
	 ExcDimensionMismatch(weights.size(), points.size()));
}



template <int dim>
Quadrature<dim>::Quadrature (const typename std::vector<Point<dim> > &points)
		:
		n_quadrature_points(points.size()),
		quadrature_points(points),
		weights(points.size(), std::atof("Inf"))
{
  Assert(weights.size() == points.size(),
	 ExcDimensionMismatch(weights.size(), points.size()));
}



template <int dim>
Quadrature<dim>::Quadrature (const Point<dim> &point):
		n_quadrature_points(1),
		quadrature_points(std::vector<Point<dim> > (1, point)),
		weights(std::vector<double> (1, 1.))
{}


template <>
Quadrature<0>::Quadrature (const Quadrature<-1> &,
			   const Quadrature<1> &)
		: n_quadrature_points (0)
{
  Assert (false, ExcInternalError());
};



template <>
Quadrature<1>::Quadrature (const Quadrature<0> &,
			   const Quadrature<1> &) :
                n_quadrature_points (0)
{
  Assert (false, ExcInternalError());
};



template <int dim>
Quadrature<dim>::Quadrature (const SubQuadrature &q1,
			     const Quadrature<1> &q2) :
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
	    = q1.point(i)(d);
	quadrature_points[present_index](dim-1)
	  = q2.point(j)(0);
					       
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
Quadrature<dim>::~Quadrature ()
{};



template <>
const Point<0> & Quadrature<0>::point (const unsigned int) const
{
  Assert (false, ExcInternalError());
  static const Point<0> dummy;
  return dummy;
};



template <int dim>
const Point<dim> & Quadrature<dim>::point (const unsigned int i) const
{
  Assert (i<n_quadrature_points, ExcIndexRange(i, 0, n_quadrature_points));
  return quadrature_points[i];
};



template <>
const std::vector<Point<0> > & Quadrature<0>::get_points () const
{
  Assert (false, ExcInternalError());
  return quadrature_points;
};



template <int dim>
const typename std::vector<Point<dim> > & Quadrature<dim>::get_points () const
{
  return quadrature_points;
};



template <>
double Quadrature<0>::weight (const unsigned int) const
{
  Assert (false, ExcInternalError());
  return 0;
};



template <int dim>
double Quadrature<dim>::weight (const unsigned int i) const
{
  Assert (i<n_quadrature_points, ExcIndexRange(i, 0, n_quadrature_points));
  return weights[i];
};



template <int dim>
const std::vector<double> & Quadrature<dim>::get_weights () const
{
  return weights;
};



template <>
const std::vector<double> & Quadrature<0>::get_weights () const
{
  Assert (false, ExcInternalError());
  return weights;
};



template <int dim>
unsigned int
Quadrature<dim>::memory_consumption () const
{
  return (MemoryConsumption::memory_consumption (quadrature_points) +
	  MemoryConsumption::memory_consumption (weights));
};


//----------------------------------------------------------------------//


template <>
void QProjector<1>::project_to_face (const Quadrature<0> &,
				     const unsigned int,
				     std::vector<Point<1> > &)
{
  Assert(false, ExcNotImplemented());
}



template <>
void QProjector<2>::project_to_face (const Quadrature<1>      &quadrature,
				     const unsigned int        face_no,
				     std::vector<Point<2> >   &q_points)
{
  const unsigned int dim=2;
  Assert (face_no<2*dim, ExcIndexRange (face_no, 0, 2*dim));
  Assert (q_points.size() == quadrature.n_quadrature_points,
	  ExcDimensionMismatch (q_points.size(), quadrature.n_quadrature_points));
  
  for (unsigned int p=0; p<quadrature.n_quadrature_points; ++p)
    switch (face_no)
      {
	case 0:
	      q_points[p] = Point<dim>(quadrature.point(p)(0),0);
	      break;	   
	case 1:
	      q_points[p] = Point<dim>(1,quadrature.point(p)(0));
	      break;	   
	case 2:
	      q_points[p] = Point<dim>(quadrature.point(p)(0),1);
	      break;	   
	case 3:
	      q_points[p] = Point<dim>(0,quadrature.point(p)(0));
	      break;
	default:
	      Assert (false, ExcInternalError());
      };
};



template <>
void QProjector<3>::project_to_face (const Quadrature<2>    &quadrature,
				     const unsigned int      face_no,
				     std::vector<Point<3> > &q_points)
{
  const unsigned int dim=3;
  Assert (face_no<2*dim, ExcIndexRange (face_no, 0, 2*dim));
  Assert (q_points.size() == quadrature.n_quadrature_points,
	  ExcDimensionMismatch (q_points.size(), quadrature.n_quadrature_points));
  
  for (unsigned int p=0; p<quadrature.n_quadrature_points; ++p)
    switch (face_no)
      {
	case 0:
	      q_points[p] = Point<dim>(quadrature.point(p)(0),
				       0,
				       quadrature.point(p)(1));
	      break;	   
	case 1:
	      q_points[p] = Point<dim>(quadrature.point(p)(0),
				       1,
				       quadrature.point(p)(1));
	      break;	   
	case 2:
	      q_points[p] = Point<dim>(quadrature.point(p)(0),
				       quadrature.point(p)(1),
				       0);
	      break;
	case 3:
	      q_points[p] = Point<dim>(1,
				       quadrature.point(p)(0),
				       quadrature.point(p)(1));
	      break;
	case 4:
	      q_points[p] = Point<dim>(quadrature.point(p)(0),
				       quadrature.point(p)(1),
				       1);
	      break;
	case 5:
	      q_points[p] = Point<dim>(0,
				       quadrature.point(p)(0),
				       quadrature.point(p)(1));
	      break;      
	      
	default:
	      Assert (false, ExcInternalError());
      };
};



template <>
void QProjector<1>::project_to_subface (const Quadrature<0> &,
					const unsigned int,
					const unsigned int,
					std::vector<Point<1> > &)
{
  Assert(false, ExcNotImplemented());
}


  
template <>
void QProjector<2>::project_to_subface (const Quadrature<1>    &quadrature,
					const unsigned int      face_no,
					const unsigned int      subface_no,
					std::vector<Point<2> > &q_points)
{
  const unsigned int dim=2;
  Assert (face_no<2*dim, ExcIndexRange (face_no, 0, 2*dim));
  Assert (subface_no<(1<<(dim-1)), ExcIndexRange (face_no, 0, 1<<(dim-1)));
  Assert (q_points.size() == quadrature.n_quadrature_points,
	  ExcDimensionMismatch (q_points.size(), quadrature.n_quadrature_points));
  
  for (unsigned int p=0; p<quadrature.n_quadrature_points; ++p)
    switch (face_no)
      {
	case 0:
	      switch (subface_no) 
		{
		  case 0:
			q_points[p]
			  = Point<dim>(quadrature.point(p)(0)/2,0);
			break;
		  case 1:
			q_points[p]
			  = Point<dim>(quadrature.point(p)(0)/2+0.5,0);
			break;
		  default:
			Assert (false, ExcInternalError());
		};
	      break;	   
	case 1:
	      switch (subface_no) 
		{
		  case 0:
			q_points[p] = Point<dim>(1,quadrature.point(p)(0)/2);
			break;
		  case 1:
			q_points[p] = Point<dim>(1,quadrature.point(p)(0)/2+0.5);
			break;
		  default:
			Assert (false, ExcInternalError());
		};
	      break;	   
	case 2:
	      switch (subface_no) 
		{
		  case 0:
			q_points[p] = Point<dim>(quadrature.point(p)(0)/2,1);
			break;
		  case 1:
			q_points[p] = Point<dim>(quadrature.point(p)(0)/2+0.5,1);
			break;
		  default:
			Assert (false, ExcInternalError());
		};
	      break;	   
	case 3:
	      switch (subface_no) 
		{
		  case 0:
			q_points[p] = Point<dim>(0,quadrature.point(p)(0)/2);
			break;
		  case 1:
			q_points[p] = Point<dim>(0,quadrature.point(p)(0)/2+0.5);
			break;
		  default:
			Assert (false, ExcInternalError());
		};
	      break;
	default:
	      Assert (false, ExcInternalError());
      };
};



template <>
void QProjector<3>::project_to_subface (const Quadrature<2>    &quadrature,
					const unsigned int      face_no,
					const unsigned int      subface_no,
					std::vector<Point<3> > &q_points)
{
  const unsigned int dim=3;
  Assert (face_no<2*dim, ExcIndexRange (face_no, 0, 2*dim));
  Assert (subface_no<(1<<(dim-1)), ExcIndexRange (face_no, 0, 1<<(dim-1)));
  Assert (q_points.size() == quadrature.n_quadrature_points,
	  ExcDimensionMismatch (q_points.size(), quadrature.n_quadrature_points));


				   // for all faces and subfaces:
				   // first project onto the first
				   // subface of each face, then move
				   // it to the right place
  for (unsigned int p=0; p<quadrature.n_quadrature_points; ++p)
    switch (face_no)
      {
	case 0:
	      q_points[p] = Point<dim>(quadrature.point(p)(0)/2,
				       0,
				       quadrature.point(p)(1)/2);
	      switch (subface_no) 
		{
		  case 0:
			break;
		  case 1:
			q_points[p][0] += 1./2.;
			break;
		  case 2:
			q_points[p][0] += 1./2.;
			q_points[p][2] += 1./2.;
			break;
		  case 3:
			q_points[p][2] += 1./2.;
			break;
		  default:
			Assert (false, ExcInternalError());
		};
	      
	      break;	   
	case 1:
	      q_points[p] = Point<dim>(quadrature.point(p)(0)/2,
				       1,
				       quadrature.point(p)(1)/2);
	      switch (subface_no) 
		{
		  case 0:
			break;
		  case 1:
			q_points[p][0] += 1./2.;
			break;
		  case 2:
			q_points[p][0] += 1./2.;
			q_points[p][2] += 1./2.;
			break;
		  case 3:
			q_points[p][2] += 1./2.;
			break;
		  default:
			Assert (false, ExcInternalError());
		};
	      break;	   
	case 2:
	      q_points[p] = Point<dim>(quadrature.point(p)(0)/2,
				       quadrature.point(p)(1)/2,
				       0);
	      switch (subface_no) 
		{
		  case 0:
			break;
		  case 1:
			q_points[p][0] += 1./2.;
			break;
		  case 2:
			q_points[p][0] += 1./2.;
			q_points[p][1] += 1./2.;
			break;
		  case 3:
			q_points[p][1] += 1./2.;
			break;
		  default:
			Assert (false, ExcInternalError());
		};
	      break;
	case 3:
	      q_points[p] = Point<dim>(1,
				       quadrature.point(p)(0)/2,
				       quadrature.point(p)(1)/2);
	      switch (subface_no) 
		{
		  case 0:
			break;
		  case 1:
			q_points[p][1] += 1./2.;
			break;
		  case 2:
			q_points[p][1] += 1./2.;
			q_points[p][2] += 1./2.;
			break;
		  case 3:
			q_points[p][2] += 1./2.;
			break;
		  default:
			Assert (false, ExcInternalError());
		};
	      break;
	case 4:
	      q_points[p] = Point<dim>(quadrature.point(p)(0)/2,
				       quadrature.point(p)(1)/2,
				       1);
	      switch (subface_no) 
		{
		  case 0:
			break;
		  case 1:
			q_points[p][0] += 1./2.;
			break;
		  case 2:
			q_points[p][0] += 1./2.;
			q_points[p][1] += 1./2.;
			break;
		  case 3:
			q_points[p][1] += 1./2.;
			break;
		  default:
			Assert (false, ExcInternalError());
		};
	      break;
	case 5:
	      q_points[p] = Point<dim>(0,
				       quadrature.point(p)(0)/2,
				       quadrature.point(p)(1)/2);
	      switch (subface_no) 
		{
		  case 0:
			break;
		  case 1:
			q_points[p][1] += 1./2.;
			break;
		  case 2:
			q_points[p][1] += 1./2.;
			q_points[p][2] += 1./2.;
			break;
		  case 3:
			q_points[p][2] += 1./2.;
			break;
		  default:
			Assert (false, ExcInternalError());
		};
	      break;      
	default:
	      Assert (false, ExcInternalError());
      };
};



template <>
Quadrature<1>
QProjector<1>::project_to_all_faces (const Quadrature<0> &)
{
  Assert (false, ExcImpossibleInDim(1));
  return Quadrature<1>(0);
};



template <int dim>
Quadrature<dim>
QProjector<dim>::project_to_all_faces (const SubQuadrature &quadrature)
{
  const unsigned int n_points = quadrature.n_quadrature_points,
		     n_faces  = GeometryInfo<dim>::faces_per_cell;

				   // first fix quadrature points
  typename std::vector<Point<dim> > q_points (n_points * n_faces);
  std::vector <Point<dim> > help(n_points);
  
				   // project to each face and copy
				   // results
  for (unsigned int face=0; face<n_faces; ++face)
    {
      project_to_face(quadrature, face, help);
      std::copy (help.begin(), help.end(), q_points.begin()+n_points*face);
    }

				   // next copy over weights
  typename std::vector<double> weights (n_points * n_faces);
  for (unsigned int face=0; face<n_faces; ++face)
    std::copy (quadrature.get_weights().begin(),
	       quadrature.get_weights().end(),
	       weights.begin()+n_points*face);
  
  return Quadrature<dim>(q_points, weights);
}



template <>
Quadrature<1>
QProjector<1>::project_to_all_subfaces (const Quadrature<0> &)
{
  Assert (false, ExcImpossibleInDim(1));
  return Quadrature<1>(0);
};



template <int dim>
Quadrature<dim>
QProjector<dim>::project_to_all_subfaces (const SubQuadrature &quadrature)
{
  const unsigned int n_points          = quadrature.n_quadrature_points,
		     n_faces           = GeometryInfo<dim>::faces_per_cell,
		     subfaces_per_face = GeometryInfo<dim>::subfaces_per_face;
  
				   // first fix quadrature points
  typename std::vector<Point<dim> > q_points (n_points * n_faces * subfaces_per_face);
  std::vector <Point<dim> > help(n_points);
  
				   // project to each face and copy
				   // results
  for (unsigned int face=0; face<n_faces; ++face)
    for (unsigned int subface=0; subface<subfaces_per_face; ++subface)
      {
	project_to_subface(quadrature, face, subface, help);
	std::copy (help.begin(), help.end(),
		   q_points.begin()+n_points*(face*subfaces_per_face+subface));
      };

				   // next copy over weights
  typename std::vector<double> weights (n_points * n_faces * subfaces_per_face);
  for (unsigned int face=0; face<n_faces; ++face)
    for (unsigned int subface=0; subface<subfaces_per_face; ++subface)
      std::copy (quadrature.get_weights().begin(),
		 quadrature.get_weights().end(),
		 weights.begin()+n_points*(face*subfaces_per_face+subface));
  
  return Quadrature<dim>(q_points, weights);
}



template <int dim>
Quadrature<dim>
QProjector<dim>::project_to_child (const Quadrature<dim>    &quadrature,
				   const unsigned int        child_no)
{
  Assert (child_no < GeometryInfo<dim>::children_per_cell,
	  ExcIndexRange (child_no, 0, GeometryInfo<dim>::children_per_cell));
  
  const unsigned int n_q_points = quadrature.n_quadrature_points;  
    
				   // projection is simple: copy over
				   // old points, scale them, and if
				   // necessary shift them
  std::vector<Point<dim> > q_points = quadrature.get_points ();
  for (unsigned int i=0; i<n_q_points; ++i)
    q_points[i] /= 2;

  switch (dim)
    {
      case 1:
      {
	if (child_no == 1)
	  for (unsigned int i=0; i<n_q_points; ++i)
	    q_points[i][0] += 0.5;
	break;
      };

      case 2:
      {
	if ((child_no == 1) || (child_no == 2))
	  for (unsigned int i=0; i<n_q_points; ++i)
	    q_points[i][0] += 0.5;
	if ((child_no == 2) || (child_no == 3))
	  for (unsigned int i=0; i<n_q_points; ++i)
	    q_points[i][1] += 0.5;
	break;
      };

      default:
	    Assert (false, ExcNotImplemented());
    };

				   // for the weights, things are
				   // equally simple: copy them and
				   // scale them
  std::vector<double> weights = quadrature.get_weights ();
  for (unsigned int i=0; i<n_q_points; ++i)
    weights[i] *= (1./GeometryInfo<dim>::children_per_cell);

  return Quadrature<dim> (q_points, weights);
};



// ------------------------------------------------------------ //


template <>
bool
QIterated<1>::uses_both_endpoints (const Quadrature<1> &base_quadrature)
{
  bool at_left = false,
      at_right = false;
  for (unsigned int i=0; i<base_quadrature.n_quadrature_points; ++i)
    {
      if (base_quadrature.point(i) == Point<1>(0.0))
	at_left = true;
      if (base_quadrature.point(i) == Point<1>(1.0))
	at_right = true;
    };

  return (at_left && at_right);
};



template <>
QIterated<1>::QIterated (const Quadrature<1> &base_quadrature,
			 const unsigned int   n_copies) :
		Quadrature<1> (uses_both_endpoints(base_quadrature) ?
			       (base_quadrature.n_quadrature_points-1) * n_copies + 1 :
			       base_quadrature.n_quadrature_points * n_copies) 
{
  Assert (n_copies >= 1, ExcInvalidNumberOfCopies(n_copies));
  
  if (!uses_both_endpoints(base_quadrature))
				     // we don't have to skip some
				     // points in order to get a
				     // reasonable quadrature formula
    {
      unsigned int next_point = 0;
      for (unsigned int copy=0; copy<n_copies; ++copy)
	for (unsigned int q_point=0; q_point<base_quadrature.n_quadrature_points; ++q_point)
	  {
	    quadrature_points[next_point] = Point<1>(base_quadrature.point(q_point)(0) / n_copies
						     +
						     (1.0*copy)/n_copies);
	    weights[next_point]           = base_quadrature.weight(q_point) / n_copies;

	    ++next_point;
	  };
    }
  else
				     // skip doubly available points
    {
      unsigned int next_point = 0;

				       // first find out the weights of
				       // the left and the right boundary
				       // points. note that these usually
				       // are but need not necessarily be
				       // the same
      double double_point_weight = 0;
      unsigned int n_end_points = 0;
      for (unsigned int i=0; i<base_quadrature.n_quadrature_points; ++i)
					 // add up the weight if this
					 // is an endpoint
	if ((base_quadrature.point(i) == Point<1>(0.0)) ||
	    (base_quadrature.point(i) == Point<1>(1.0)))
	  {
	    double_point_weight += base_quadrature.weight(i);
	    ++n_end_points;
	  };
				       // scale the weight correctly
      double_point_weight /= n_copies;

				       // make sure the base quadrature formula
				       // has only one quadrature point
				       // per end point
      Assert (n_end_points == 2, ExcInvalidQuadratureFormula());


      for (unsigned int copy=0; copy<n_copies; ++copy)
	for (unsigned int q_point=0; q_point<base_quadrature.n_quadrature_points; ++q_point)
	  {
					     // skip the left point of
					     // this copy since we
					     // have already entered
					     // it the last time
	    if ((copy > 0) &&
		(base_quadrature.point(q_point) == Point<1>(0.0)))
	      continue;
	    
	    quadrature_points[next_point] = Point<1>(base_quadrature.point(q_point)(0) / n_copies
						     +
						     (1.0*copy)/n_copies);

					     // if this is the
					     // rightmost point of one
					     // of the non-last
					     // copies: give it the
					     // double weight
	    if ((copy != n_copies-1) &&
		(base_quadrature.point(q_point) == Point<1>(1.0)))
	      weights[next_point] = double_point_weight;
	    else
	      weights[next_point] = base_quadrature.weight(q_point) / n_copies;
	    
	    ++next_point;
	  };
    };

#if DEBUG
  double sum_of_weights = 0;
  for (unsigned int i=0; i<n_quadrature_points; ++i)
    sum_of_weights += weight(i);
  Assert (std::fabs(sum_of_weights-1) < 1e-15,
	  ExcSumOfWeightsNotOne());
#endif
};



// construct higher dimensional quadrature formula by tensor product
// of lower dimensional iterated quadrature formulae
template <int dim>
QIterated<dim>::QIterated (const Quadrature<1> &base_quadrature,
			   const unsigned int   N) :
		Quadrature<dim> (QIterated<dim-1>(base_quadrature, N),
				 QIterated<1>(base_quadrature, N))
{};



// explicit instantiations; note: we need them all for all dimensions
template class Quadrature<1>;
template class Quadrature<2>;
template class Quadrature<3>;
template class QIterated<1>;
template class QIterated<2>;
template class QIterated<3>;
template class QProjector<1>;
template class QProjector<2>;
template class QProjector<3>;
