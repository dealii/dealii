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

//TODO: Discussion on an efficient implementation of Point additions.

template <int dim>
double
FunctionDerivative<dim>::value (const Point<dim>   &p,
				const unsigned int  component) const
{
  switch (formula)
    {
      case Euler:
	    return (f.value(p+incr, component)-f.value(p-incr, component))/(2*h);
      case UpwindEuler:
	    return (f.value(p, component)-f.value(p-incr, component))/h;
      case FourthOrder:
	    return (-f.value(p+2*incr, component) + 8*f.value(p+incr, component)
		    -8*f.value(p-incr, component) + f.value(p-2*incr, component))/(12*h);
      default:
	    Assert(false, ExcInvalidFormula());
    }
  return 0.;
}



template <int dim>
void
FunctionDerivative<dim>::value_list (const std::vector<Point<dim> > &points,
				     std::vector<double>            &values,
				     const unsigned int              component) const
{
  const unsigned int n = points.size();
  
  switch (formula)
    {
      case Euler:
      {
	std::vector<Point<dim> > p1(n);
	std::vector<Point<dim> > p2(n);
	std::vector<double> e2(n);
	for (unsigned int i=0;i<n;++i)
	  {
	    p1[i] = points[i]+incr;
	    p2[i] = points[i]-incr;
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
	  p2[i] = points[i]-incr;
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
	    p_p[i] = points[i]+incr;
	    p_pp[i] = p_p[i]+incr;
	    p_m[i] = points[i]-incr;
	    p_mm[i] = p_m[i]-incr;
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


template FunctionDerivative<1>;
template FunctionDerivative<2>;
template FunctionDerivative<3>;
