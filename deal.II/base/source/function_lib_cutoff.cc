//----------------------------  function_lib.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002 by the deal authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  function_lib.cc  ---------------------------


#include <base/tensor.h>
#include <base/point.h>
#include <base/function_lib.h>

#include <cmath>


// in strict ANSI C mode, the following constants are not defined by
// default, so we do it ourselves
//TODO:[?] Unify the various places where PI is defined to a central instance 
#ifndef M_PI
#  define	M_PI		3.14159265358979323846
#endif

#ifndef M_PI_2
#  define	M_PI_2		1.57079632679489661923
#endif

#ifndef M_E
#  define       M_E             2.7182818284590452354
#endif

namespace Functions
{

  using namespace std;
  
  template<int dim>
  CutOffFunctionLinfty<dim>::CutOffFunctionLinfty (const double r,
						   const Point<dim> p)
		  :
		  Function<dim> (1),
                  center(p),
                  radius(r)
  {}


  template<int dim>
  double
  CutOffFunctionLinfty<dim>::value (const Point<dim>   &p,
				    const unsigned int) const
  {
    return ((center.distance(p)<radius) ? 1. : 0.);
  }


  template<int dim>
  void 
  CutOffFunctionLinfty<dim>::value_list (const typename std::vector<Point<dim> > &points,
					 std::vector<double>                     &values,
					 const unsigned int) const
  {
    Assert (values.size() == points.size(),
	    ExcDimensionMismatch(values.size(), points.size()));

    for (unsigned int i=0;i<values.size();++i)
      values[i] = (center.distance(points[i])<radius) ? 1. : 0.;
  }



  template<int dim>
  CutOffFunctionW1<dim>::CutOffFunctionW1 (const double     r,
					   const Point<dim> p)
		  :
		  Function<dim> (1),
                  center(p),
                  radius(r)
  {}


  template<int dim>
  double
  CutOffFunctionW1<dim>::value (const Point<dim>   &p,
				const unsigned int) const
  {
    const double d = center.distance(p);
    return ((d<radius) ? (radius-d) : 0.);
  }


  template<int dim>
  void 
  CutOffFunctionW1<dim>::value_list (const typename std::vector<Point<dim> > &points,
				     std::vector<double>                     &values,
				     const unsigned int) const
  {
    Assert (values.size() == points.size(),
	    ExcDimensionMismatch(values.size(), points.size()));
    
    for (unsigned int i=0;i<values.size();++i)
      {
	const double d = center.distance(points[i]);
	values[i] = ((d<radius) ? (radius-d) : 0.);
      }
  }



  template<int dim>
  CutOffFunctionCinfty<dim>::CutOffFunctionCinfty (const double     r,
						   const Point<dim> p)
		  :
		  Function<dim> (1),
                  center(p),
                  radius(r)
  {}


  template<int dim>
  double
  CutOffFunctionCinfty<dim>::value (const Point<dim>   &p,
				    const unsigned int) const
  {
    const double d = center.distance(p);
    const double r = radius;
    if (d>=r)
      return 0.;
    const double e = -r*r/(r*r-d*d);
    return  ((e<-50) ? 0. : M_E * exp(e));
  }


  template<int dim>
  void 
  CutOffFunctionCinfty<dim>::value_list (const typename std::vector<Point<dim> > &points,
					 std::vector<double>                     &values,
					 const unsigned int) const
  {
    Assert (values.size() == points.size(),
	    ExcDimensionMismatch(values.size(), points.size()));
    
    const double r = radius;

    for (unsigned int i=0;i<values.size();++i)
      {
	const double d = center.distance(points[i]);
	if (d>=r)
	  {
	    values[i] = 0.;
	  } else {
	    const double e = -r*r/(r*r-d*d);
	    values[i] = (e<-50) ? 0. : M_E * exp(e);
	  }
      }
  }



  template<int dim>
  Tensor<1,dim>
  CutOffFunctionCinfty<dim>::gradient (const Point<dim>   &p,
				       const unsigned int) const
  {
    const double d = center.distance(p);
    const double r = radius;
    if (d>=r)
      return Tensor<1,dim>();
    const double e = -d*d/(r-d)/(r+d);
    return  ((e<-50) ?
	     Point<dim>() :
	     (p-center)/d*(-2.0*r*r/pow(-r*r+d*d,2.0)*d*exp(e)));
  }
  

// explicit instantiations  
  template class CutOffFunctionLinfty <1>;
  template class CutOffFunctionLinfty <2>;
  template class CutOffFunctionLinfty <3>;
  
  template class CutOffFunctionW1 <1>;
  template class CutOffFunctionW1 <2>;
  template class CutOffFunctionW1 <3>;
  
  template class CutOffFunctionCinfty <1>;
  template class CutOffFunctionCinfty <2>;
  template class CutOffFunctionCinfty <3>;
}
