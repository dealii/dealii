// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------


#include <deal.II/base/geometry_info.h>
#include <deal.II/base/polynomials_rannacher_turek.h>

#include <memory>

DEAL_II_NAMESPACE_OPEN


template <int dim>
PolynomialsRannacherTurek<dim>::PolynomialsRannacherTurek()
  : ScalarPolynomialsBase<dim>(2, dealii::GeometryInfo<dim>::faces_per_cell)
{
  Assert(dim == 2, ExcNotImplemented());
}



template <int dim>
double
PolynomialsRannacherTurek<dim>::compute_value(const unsigned int i,
                                              const Point<dim> & p) const
{
  Assert(dim == 2, ExcNotImplemented());
  if (i == 0)
    {
      return (0.75 - 2.5 * p(0) + 1.5 * p(1) +
              1.5 * (p(0) * p(0) - p(1) * p(1)));
    }
  else if (i == 1)
    {
      return (-0.25 - 0.5 * p(0) + 1.5 * p(1) +
              1.5 * (p(0) * p(0) - p(1) * p(1)));
    }
  else if (i == 2)
    {
      return (0.75 + 1.5 * p(0) - 2.5 * p(1) -
              1.5 * (p(0) * p(0) - p(1) * p(1)));
    }
  else if (i == 3)
    {
      return (-0.25 + 1.5 * p(0) - 0.5 * p(1) -
              1.5 * (p(0) * p(0) - p(1) * p(1)));
    }

  Assert(false, ExcNotImplemented());
  return 0;
}



template <int dim>
Tensor<1, dim>
PolynomialsRannacherTurek<dim>::compute_grad(const unsigned int i,
                                             const Point<dim> & p) const
{
  Assert(dim == 2, ExcNotImplemented());
  Tensor<1, dim> grad;
  if (i == 0)
    {
      grad[0] = -2.5 + 3 * p(0);
      grad[1] = 1.5 - 3 * p(1);
    }
  else if (i == 1)
    {
      grad[0] = -0.5 + 3.0 * p(0);
      grad[1] = 1.5 - 3.0 * p(1);
    }
  else if (i == 2)
    {
      grad[0] = 1.5 - 3.0 * p(0);
      grad[1] = -2.5 + 3.0 * p(1);
    }
  else if (i == 3)
    {
      grad[0] = 1.5 - 3.0 * p(0);
      grad[1] = -0.5 + 3.0 * p(1);
    }
  else
    {
      Assert(false, ExcNotImplemented());
    }

  return grad;
}



template <int dim>
Tensor<2, dim>
PolynomialsRannacherTurek<dim>::compute_grad_grad(
  const unsigned int i,
  const Point<dim> & /*p*/) const
{
  Assert(dim == 2, ExcNotImplemented());
  Tensor<2, dim> grad_grad;
  if (i == 0)
    {
      grad_grad[0][0] = 3;
      grad_grad[0][1] = 0;
      grad_grad[1][0] = 0;
      grad_grad[1][1] = -3;
    }
  else if (i == 1)
    {
      grad_grad[0][0] = 3;
      grad_grad[0][1] = 0;
      grad_grad[1][0] = 0;
      grad_grad[1][1] = -3;
    }
  else if (i == 2)
    {
      grad_grad[0][0] = -3;
      grad_grad[0][1] = 0;
      grad_grad[1][0] = 0;
      grad_grad[1][1] = 3;
    }
  else if (i == 3)
    {
      grad_grad[0][0] = -3;
      grad_grad[0][1] = 0;
      grad_grad[1][0] = 0;
      grad_grad[1][1] = 3;
    }
  return grad_grad;
}



template <int dim>
void
PolynomialsRannacherTurek<dim>::evaluate(
  const Point<dim> &           unit_point,
  std::vector<double> &        values,
  std::vector<Tensor<1, dim>> &grads,
  std::vector<Tensor<2, dim>> &grad_grads,
  std::vector<Tensor<3, dim>> &third_derivatives,
  std::vector<Tensor<4, dim>> &fourth_derivatives) const
{
  const unsigned int n_pols = this->n();
  Assert(values.size() == n_pols || values.size() == 0,
         ExcDimensionMismatch(values.size(), n_pols));
  Assert(grads.size() == n_pols || grads.size() == 0,
         ExcDimensionMismatch(grads.size(), n_pols));
  Assert(grad_grads.size() == n_pols || grad_grads.size() == 0,
         ExcDimensionMismatch(grad_grads.size(), n_pols));
  Assert(third_derivatives.size() == n_pols || third_derivatives.size() == 0,
         ExcDimensionMismatch(third_derivatives.size(), n_pols));
  Assert(fourth_derivatives.size() == n_pols || fourth_derivatives.size() == 0,
         ExcDimensionMismatch(fourth_derivatives.size(), n_pols));

  for (unsigned int i = 0; i < n_pols; ++i)
    {
      if (values.size() != 0)
        {
          values[i] = compute_value(i, unit_point);
        }
      if (grads.size() != 0)
        {
          grads[i] = compute_grad(i, unit_point);
        }
      if (grad_grads.size() != 0)
        {
          grad_grads[i] = compute_grad_grad(i, unit_point);
        }
      if (third_derivatives.size() != 0)
        {
          third_derivatives[i] = compute_derivative<3>(i, unit_point);
        }
      if (fourth_derivatives.size() != 0)
        {
          fourth_derivatives[i] = compute_derivative<4>(i, unit_point);
        }
    }
}



template <int dim>
std::unique_ptr<ScalarPolynomialsBase<dim>>
PolynomialsRannacherTurek<dim>::clone() const
{
  return std::make_unique<PolynomialsRannacherTurek<dim>>(*this);
}


// explicit instantiations
#include "polynomials_rannacher_turek.inst"

DEAL_II_NAMESPACE_CLOSE
