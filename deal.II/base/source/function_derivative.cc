//--------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000 by the deal.II authors
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
FunctionDerivative<dim>::FunctionDerivative (const Function<dim>& f,
					     const Point<dim>& dir)
  :
  Function<dim> (f.n_components, f.get_time()),
  f(f),
  direction(dir)
{
  set_h();
  set_formula();
}



template <int dim>
void
FunctionDerivative<dim>::set_h (double h_new)
{
  h = h_new;
  incr = h*direction;
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
  switch (formula)
    {
    case Euler:
      return (f.value(p+incr, component)-f.value(p-incr, component))/2*h;
    case UpwindEuler:
      return (f.value(p, component)-f.value(p-incr, component))/h;
    default:
      Assert(false, ExcInvalidFormula());
    }
  return 0.;
}



template <int dim>
void
FunctionDerivative<dim>::value_list (const vector<Point<dim> > &points,
				     vector<double>            &values,
				     const unsigned int         component) const
{
  const unsigned int n = points.size();
  
  switch (formula)
    {
    case Euler:
    {
      vector<Point<dim> > p1(n);
      vector<Point<dim> > p2(n);
      vector<double> e2(n);
      for (unsigned int i=0;i<n;++i)
	{
	  p1[i] = points[i]+incr;
	  p2[i] = points[i]-incr;
	}
      f.value_list(p1, values, component);
      f.value_list(p2, e2, component);
      for (unsigned int i=0;i<n;++i)
	values[i] = (values[i]-e2[i])/2*h;
      break;
    }    
    case UpwindEuler:
    {
      vector<Point<dim> > p2(n);
      vector<double> e2(n);
      for (unsigned int i=0;i<n;++i)
	p2[i] = points[i]-incr;
      f.value_list(points, values, component);
      f.value_list(p2, e2, component);
      for (unsigned int i=0;i<n;++i)
	values[i] = (values[i]-e2[i])/h;
      break;
    }
    default:
      Assert(false, ExcInvalidFormula());
    }
}


template FunctionDerivative<1>;
template FunctionDerivative<2>;
template FunctionDerivative<3>;
