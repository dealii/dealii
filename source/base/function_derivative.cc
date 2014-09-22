// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2013 by the deal.II authors
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

#include <deal.II/base/point.h>
#include <deal.II/base/function_derivative.h>
#include <deal.II/lac/vector.h>

#include <cmath>

DEAL_II_NAMESPACE_OPEN

template <int dim>
FunctionDerivative<dim>::FunctionDerivative (const Function<dim> &f,
                                             const Point<dim>    &dir,
                                             const double         h)
  :
  AutoDerivativeFunction<dim> (h, f.n_components, f.get_time()),
  f(f),
  h(h),
  incr(1, h *dir)
{
  set_formula();
}



template <int dim>
FunctionDerivative<dim>::FunctionDerivative (const Function<dim> &f,
                                             const std::vector<Point<dim> > &dir,
                                             const double h)
  :
  AutoDerivativeFunction<dim> (h, f.n_components, f.get_time()),
  f(f),
  h(h),
  incr(dir.size())
{
  for (unsigned int i=0; i<incr.size (); ++i)
    incr[i] = h*dir[i];
  set_formula();
}



template <int dim>
void
FunctionDerivative<dim>::set_formula (typename AutoDerivativeFunction<dim>::DifferenceFormula form)
{
  formula = form;
}



template <int dim>
void
FunctionDerivative<dim>::set_h (const double new_h)
{
  for (unsigned int i=0; i<incr.size (); ++i)
    incr[i] *= new_h/h;
  h = new_h;
}



template <int dim>
double
FunctionDerivative<dim>::value (const Point<dim>   &p,
                                const unsigned int  component) const
{
  Assert (incr.size() == 1,
          ExcMessage ("FunctionDerivative was not initialized for constant direction"));

  switch (formula)
    {
    case AutoDerivativeFunction<dim>::Euler:
      return (f.value(p+incr[0], component)-f.value(p-incr[0], component))/(2*h);
    case AutoDerivativeFunction<dim>::UpwindEuler:
      return (f.value(p, component)-f.value(p-incr[0], component))/h;
    case AutoDerivativeFunction<dim>::FourthOrder:
      return (-f.value(p+2*incr[0], component) + 8*f.value(p+incr[0], component)
              -8*f.value(p-incr[0], component) + f.value(p-2*incr[0], component))/(12*h);
    default:
      Assert(false, ExcInvalidFormula());
    }
  return 0.;
}



template <int dim>
void
FunctionDerivative<dim>::vector_value (
  const Point<dim>   &p,
  Vector<double> &result) const
{
  Assert (incr.size() == 1,
          ExcMessage ("FunctionDerivative was not initialized for constant direction"));
  Vector<double> aux(result.size());

  // Formulas are the same as in
  // value, but here we have to use
  // Vector arithmetic
  switch (formula)
    {
    case AutoDerivativeFunction<dim>::Euler:
      f.vector_value(p+incr[0], result);
      f.vector_value(p-incr[0], aux);
      result.sadd(1./(2*h), -1./(2*h), aux);
      return;
    case AutoDerivativeFunction<dim>::UpwindEuler:
      f.vector_value(p, result);
      f.vector_value(p-incr[0], aux);
      result.sadd(1./h, -1./h, aux);
      return;
    case AutoDerivativeFunction<dim>::FourthOrder:
      f.vector_value(p-2*incr[0], result);
      f.vector_value(p+2*incr[0], aux);
      result.add(-1., aux);
      f.vector_value(p-incr[0], aux);
      result.add(-8., aux);
      f.vector_value(p+incr[0], aux);
      result.add(8., aux);
      result/=(12.*h);
      return;
    default:
      Assert(false, ExcInvalidFormula());
    }
}



template <int dim>
void
FunctionDerivative<dim>::value_list (const std::vector<Point<dim> > &points,
                                     std::vector<double>            &values,
                                     const unsigned int              component) const
{
  const unsigned int n = points.size();
  const bool variable_direction = (incr.size() == 1) ? false : true;
  if (variable_direction)
    Assert (incr.size() == points.size(),
            ExcDimensionMismatch(incr.size(), points.size()));

  // Vector of auxiliary values
  std::vector<double> aux(n);
  // Vector of auxiliary points
  std::vector<Point<dim> > paux(n);
  // Use the same formulas as in
  // value, but with vector
  // arithmetic
  switch (formula)
    {
    case AutoDerivativeFunction<dim>::Euler:
      for (unsigned int i=0; i<n; ++i)
        paux[i] = points[i]+incr[(variable_direction) ? i : 0U];
      f.value_list(paux, values, component);
      for (unsigned int i=0; i<n; ++i)
        paux[i] = points[i]-incr[(variable_direction) ? i : 0U];
      f.value_list(paux, aux, component);
      for (unsigned int i=0; i<n; ++i)
        values[i] = (values[i]-aux[i])/(2*h);
      return;
    case AutoDerivativeFunction<dim>::UpwindEuler:
      f.value_list(points, values, component);
      for (unsigned int i=0; i<n; ++i)
        paux[i] = points[i]-incr[(variable_direction) ? i : 0U];
      f.value_list(paux, aux, component);
      for (unsigned int i=0; i<n; ++i)
        values[i] = (values[i]-aux[i])/h;
      return;
    case AutoDerivativeFunction<dim>::FourthOrder:
      for (unsigned int i=0; i<n; ++i)
        paux[i] = points[i]-2*incr[(variable_direction) ? i : 0U];
      f.value_list(paux, values, component);
      for (unsigned int i=0; i<n; ++i)
        paux[i] = points[i]+2*incr[(variable_direction) ? i : 0U];
      f.value_list(paux, aux, component);
      for (unsigned int i=0; i<n; ++i)
        values[i] -= aux[i];
      for (unsigned int i=0; i<n; ++i)
        paux[i] = points[i]+incr[(variable_direction) ? i : 0U];
      f.value_list(paux, aux, component);
      for (unsigned int i=0; i<n; ++i)
        values[i] += 8.*aux[i];
      for (unsigned int i=0; i<n; ++i)
        paux[i] = points[i]-incr[(variable_direction) ? i : 0U];
      f.value_list(paux, aux, component);
      for (unsigned int i=0; i<n; ++i)
        values[i] = (values[i] - 8.*aux[i])/(12*h);
      return;
    default:
      Assert(false, ExcInvalidFormula());
    }
}



template <int dim>
std::size_t
FunctionDerivative<dim>::memory_consumption () const
{
  // only simple data elements, so
  // use sizeof operator
  return sizeof (*this);
}



// explicit instantiations
template class FunctionDerivative<1>;
template class FunctionDerivative<2>;
template class FunctionDerivative<3>;

DEAL_II_NAMESPACE_CLOSE
