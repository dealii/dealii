//--------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001 by the deal authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//--------------------------------------------------------------------------


#include <base/point.h>
#include <base/function_derivative.h>

#include <cmath>

template <int dim>
FunctionDerivative<dim>::FunctionDerivative (const Function<dim> &f,
					     const Point<dim>    &dir,
					     const double         h)
		:
		Function<dim> (f.n_components, f.get_time()),
                f(f),
                h(h),
                incr(1, h*dir)
{
  set_formula();
}



template <int dim>
FunctionDerivative<dim>::FunctionDerivative (const Function<dim>& f,
					     const std::vector<Point<dim> >& dir,
					     const double h)
		:
		Function<dim> (f.n_components, f.get_time()),
  f(f),
  h(h),
  incr(dir.size())
{
  for (unsigned int i=0;i<incr.size ();++i)
    incr[i] = h*dir[i];
  set_formula();
}



template <int dim>
void
FunctionDerivative<dim>::set_formula (DifferenceFormula form)
{
  formula = form;
}

//TODO:[?] Discussion on an efficient implementation of Point additions.

template <int dim>
double
FunctionDerivative<dim>::value (const Point<dim>   &p,
				const unsigned int  component) const
{
  Assert (incr.size() == 1,
	  ExcMessage ("FunctionDerivative was not initialized for constant direection"));

  switch (formula)
    {
      case Euler:
	    return (f.value(p+incr[0], component)-f.value(p-incr[0], component))/(2*h);
      case UpwindEuler:
	    return (f.value(p, component)-f.value(p-incr[0], component))/h;
      case FourthOrder:
	    return (-f.value(p+2*incr[0], component) + 8*f.value(p+incr[0], component)
		    -8*f.value(p-incr[0], component) + f.value(p-2*incr[0], component))/(12*h);
      default:
	    Assert(false, ExcInvalidFormula());
    }
  return 0.;
}



// if necessary try to work around a bug in the IBM xlC compiler
#ifdef XLC_WORK_AROUND_STD_BUG
using namespace std;
#endif

//TODO:[?] Optimize construction of vectors thread-safe

template <int dim>
void
FunctionDerivative<dim>::value_list (const typename std::vector<Point<dim> > &points,
				     std::vector<double>            &values,
				     const unsigned int              component) const
{
  const unsigned int n = points.size();
  bool variable_direction = (incr.size() == 1) ? false : true;
  if (variable_direction)
    Assert (incr.size() == points.size(),
	    ExcDimensionMismatch(incr.size(), points.size()));

  switch (formula)
    {
      case Euler:
      {
	std::vector<Point<dim> > p1(n);
	std::vector<Point<dim> > p2(n);
	std::vector<double> e2(n);
	for (unsigned int i=0;i<n;++i)
	  {
	    const unsigned int j = (variable_direction) ? i:0;
	    p1[i] = points[i]+incr[j];
	    p2[i] = points[i]-incr[j];
	  }
	f.value_list(p1, values, component);
	f.value_list(p2, e2, component);
      
	for (unsigned int i=0;i<n;++i)
	  {
	    values[i] = (values[i]-e2[i])/(2*h);
	  }
	break;
      }    
      case UpwindEuler:
      {
	std::vector<Point<dim> > p2(n);
	std::vector<double> e2(n);
	for (unsigned int i=0;i<n;++i)
	  {
	    const unsigned int j = (variable_direction) ? i:0;
	    p2[i] = points[i]-incr[j];
	  }
	f.value_list(points, values, component);
	f.value_list(p2, e2, component);
	for (unsigned int i=0;i<n;++i)
	  values[i] = (values[i]-e2[i])/h;
	break;
      }
      case FourthOrder:
      {
	std::vector<Point<dim> > p_p(n);
	std::vector<Point<dim> > p_pp(n);
	std::vector<Point<dim> > p_m(n);
	std::vector<Point<dim> > p_mm(n);
	std::vector<double> e_p(n);
	std::vector<double> e_pp(n);
	std::vector<double> e_m(n);
	for (unsigned int i=0;i<n;++i)
	  {
	    const unsigned int j = (variable_direction) ? i:0;
	    p_p[i] = points[i]+incr[j];
	    p_pp[i] = p_p[i]+incr[j];
	    p_m[i] = points[i]-incr[j];
	    p_mm[i] = p_m[i]-incr[j];
	  }
	f.value_list(p_mm, values, component);
	f.value_list(p_pp, e_pp, component);
	f.value_list(p_p, e_p, component);
	f.value_list(p_m, e_m, component);
      
	for (unsigned int i=0;i<n;++i)
	  {
	    values[i] = (values[i]-e_pp[i]+8*(e_p[i]-e_m[i]))/(12*h);
	  }
	break;
      }    

      default:
	    Assert(false, ExcInvalidFormula());
    }
}



template <int dim>
unsigned int
FunctionDerivative<dim>::memory_consumption () const
{
				   // only simple data elements, so
				   // use sizeof operator
  return sizeof (*this);
};


template class FunctionDerivative<1>;
template class FunctionDerivative<2>;
template class FunctionDerivative<3>;
