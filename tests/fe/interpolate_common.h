//----------------------------------------------------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2005, 2006, 2008, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------------------------------------------------

#include "../tests.h"
#include <base/quadrature_lib.h>
#include <base/function.h>
#include <lac/vector.h>
#include <fe/fe.h>

#include <iomanip>

// Compute the maximal difference between local finite element interpolant and
// a given polynomial function. If the finite element space is rich enough,
// then the error should be zero

template <int dim>
double difference(
  const FiniteElement<dim>& fe,
  const std::vector<double> dofs,
  const Function<dim>& function)
{
  double result = 0.;
  QGauss<dim> quadrature(fe.degree+1);
  
  std::vector<double> f(quadrature.size());
  function.value_list(quadrature.get_points(), f);
  
  for (unsigned int k=0;k<quadrature.size();++k)
    {
      double diff = f[k];
      for (unsigned int i=0;i<dofs.size();++i)
	diff -= dofs[i] * fe.shape_value(i, quadrature.point(k));
      diff = std::abs(diff);
      result = std::max(result, diff);
    }
  return result;
}

template <int dim>
double vector_difference(
  const FiniteElement<dim>& fe,
  const std::vector<double> dofs,
  const Function<dim>& function,
  const unsigned int offset)
{
  double result = 0.;
  QGauss<dim> quadrature(fe.degree+1);
  
  std::vector<Vector<double> > f(quadrature.size(),
				 Vector<double>(function.n_components));

  function.vector_value_list(quadrature.get_points(), f);
  
  for (unsigned int k=0;k<quadrature.size();++k)
    for (unsigned int comp=0;comp<fe.n_components();++comp)
      {
	double diff = f[k](comp+offset);
	for (unsigned int i=0;i<dofs.size();++i)
	  diff -= dofs[i] * fe.shape_value_component(i, quadrature.point(k),
						     comp);
	diff = std::abs(diff);
	result = std::max(result, diff);
      }
  return result;
}


// Local implementation for any dimension


template<int dim, int degree, int COMP=1>
class
Q1WedgeFunction : public Function<dim>
{
  public:
    Q1WedgeFunction() : Function<dim> (COMP)
      {}
    
    double value (const Point<dim>   &p,
		  const unsigned int c) const
      {
	double result = 1.;
	for (unsigned int d=0;d<dim;++d)
	  for (unsigned int k=0;k<degree;++k)
	    result *= p(d) + c;
	return result;
      }
    
    void value_list (const std::vector<Point<dim> > &points,
		     std::vector<double>            &values,
		     const unsigned int c) const
      {
	Assert (values.size() == points.size(),
		ExcDimensionMismatch(values.size(), points.size()));
	
	for (unsigned int i=0;i<points.size();++i)
	  {
	    const Point<dim>& p = points[i];
	    double result = 1.;
	    for (unsigned int d=0;d<dim;++d)
	      for (unsigned int k=0;k<degree;++k)
		result *= p(d) + c;
	    values[i] = result;
	  }
      }

    void vector_value_list (const std::vector<Point<dim> > &points,
			    std::vector<Vector<double> >&   values) const
      {
	Assert (values.size() == points.size(),
		ExcDimensionMismatch(values.size(), points.size()));
	Assert (values[0].size() == this->n_components,
		ExcDimensionMismatch(values.size(), this->n_components));
	
	for (unsigned int i=0;i<points.size();++i)
	  {
	    const Point<dim>& p = points[i];
	    for (unsigned int c=0;c<COMP;++c)
	      {
		double result = 1.;
		for (unsigned int d=0;d<dim;++d)
		  for (unsigned int k=0;k<degree;++k)
		    result *= p(d);
		values[i](c) = result;
	      }
	  }
      }
};

