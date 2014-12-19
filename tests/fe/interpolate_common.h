// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


#include "../tests.h"
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/lac/vector.h>
#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_q_dg0.h>

#include <iomanip>

// Compute the maximal difference between local finite element interpolant and
// a given polynomial function. If the finite element space is rich enough,
// then the error should be zero

template <int dim>
double difference(
  const FiniteElement<dim> &fe,
  const std::vector<double> dofs,
  const Function<dim> &function)
{
  double result = 0.;
  QGauss<dim> quadrature(fe.degree+1);

  std::vector<double> f(quadrature.size());
  function.value_list(quadrature.get_points(), f);

  for (unsigned int k=0; k<quadrature.size(); ++k)
    {
      double diff = f[k];
      for (unsigned int i=0; i<dofs.size(); ++i)
        diff -= dofs[i] * fe.shape_value(i, quadrature.point(k));
      diff = std::abs(diff);
      result = std::max(result, diff);
    }
  return result;
}

template <int dim>
double vector_difference(
  const FiniteElement<dim> &fe,
  const std::vector<double> dofs,
  const Function<dim> &function,
  const unsigned int offset)
{
  double result = 0.;
  QGauss<dim> quadrature(fe.degree+1);

  std::vector<Vector<double> > f(quadrature.size(),
                                 Vector<double>(function.n_components));

  function.vector_value_list(quadrature.get_points(), f);

  for (unsigned int k=0; k<quadrature.size(); ++k)
    for (unsigned int comp=0; comp<fe.n_components(); ++comp)
      {
        double diff = f[k](comp+offset);
        for (unsigned int i=0; i<dofs.size(); ++i)
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
    for (unsigned int d=0; d<dim; ++d)
      for (unsigned int k=0; k<degree; ++k)
        result *= p(d) + c;
    return result;
  }

  void value_list (const std::vector<Point<dim> > &points,
                   std::vector<double>            &values,
                   const unsigned int c) const
  {
    Assert (values.size() == points.size(),
            ExcDimensionMismatch(values.size(), points.size()));

    for (unsigned int i=0; i<points.size(); ++i)
      {
        const Point<dim> &p = points[i];
        double result = 1.;
        for (unsigned int d=0; d<dim; ++d)
          for (unsigned int k=0; k<degree; ++k)
            result *= p(d) + c;
        values[i] = result;
      }
  }

  void vector_value_list (const std::vector<Point<dim> > &points,
                          std::vector<Vector<double> >   &values) const
  {
    Assert (values.size() == points.size(),
            ExcDimensionMismatch(values.size(), points.size()));
    Assert (values[0].size() == this->n_components,
            ExcDimensionMismatch(values.size(), this->n_components));

    for (unsigned int i=0; i<points.size(); ++i)
      {
        const Point<dim> &p = points[i];
        for (unsigned int c=0; c<COMP; ++c)
          {
            double result = 1.;
            for (unsigned int d=0; d<dim; ++d)
              for (unsigned int k=0; k<degree; ++k)
                result *= p(d);
            values[i](c) = result;
          }
      }
  }
};

