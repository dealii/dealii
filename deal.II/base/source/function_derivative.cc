//--------------------------------------------------------------------------
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
					     const typename std::vector<Point<dim> >& dir,
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



template <int dim>
void
FunctionDerivative<dim>::value_list (const typename std::vector<Point<dim> > &points,
				     std::vector<double>            &values,
				     const unsigned int              component) const
{
  const unsigned int n = points.size();
  const bool variable_direction = (incr.size() == 1) ? false : true;
  if (variable_direction)
    Assert (incr.size() == points.size(),
	    ExcDimensionMismatch(incr.size(), points.size()));

  switch (formula)
    {
      case Euler:
      {
					 // let p1 and p2 be arrays of
					 // evaluation points shifted
					 // a little in direction j
	std::vector<Point<dim> > p1 = points;
	std::vector<Point<dim> > p2 = points;

	for (unsigned int i=0; i<n; ++i)
	  {
	    const unsigned int j = (variable_direction) ? i:0;
	    p1[i] += incr[j];
	    p2[i] -= incr[j];
	  }	

					 // next get values of
					 // functions at these points
	std::vector<double> values2(n);
	f.value_list(p1, values,  component);
	f.value_list(p2, values2, component);

					 // finally compute finite
					 // differences
	for (unsigned int i=0; i<n; ++i)
	  values[i] = (values[i]-values2[i])/(2*h);
	  
	break;
      }
       
      case UpwindEuler:
      {
					 // compute upwind points
	std::vector<Point<dim> > p2 = points;
	for (unsigned int i=0; i<n; ++i)
	  {
	    const unsigned int j = (variable_direction) ? i:0;
	    p2[i] -= incr[j];
	  }

					 // get values at points
	std::vector<double> values2(n);
	f.value_list(points, values,  component);
	f.value_list(p2,     values2, component);

					 // compute finite differences
	for (unsigned int i=0; i<n; ++i)
	  values[i] = (values[i]-values2[i])/h;
	break;
      }
       
      case FourthOrder:
      {
					 // first compute evaluation
					 // points
	std::vector<Point<dim> > p_p = points;
	std::vector<Point<dim> > p_pp(n);
	std::vector<Point<dim> > p_m = points;
	std::vector<Point<dim> > p_mm(n);
	for (unsigned int i=0;i<n;++i)
	  {
	    const unsigned int j = (variable_direction) ? i:0;
	    p_p[i] += incr[j];
	    p_pp[i] = p_p[i]+incr[j];
	    p_m[i] -= incr[j];
	    p_mm[i] = p_m[i]-incr[j];
	  }

					 // next compute values of
					 // function at these
					 // points. use @p{values} for
					 // @p{e_mm}
	std::vector<double> e_p(n);
	std::vector<double> e_pp(n);
	std::vector<double> e_m(n);

	f.value_list(p_mm, values, component);
	f.value_list(p_pp, e_pp,   component);
	f.value_list(p_p,  e_p,    component);
	f.value_list(p_m,  e_m,    component);

					 // compute finite differences
	for (unsigned int i=0; i<n; ++i)
	  values[i] = (values[i]-e_pp[i]+8*(e_p[i]-e_m[i]))/(12*h);
	  
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



// explicit instantiations
template class FunctionDerivative<1>;
template class FunctionDerivative<2>;
template class FunctionDerivative<3>;
