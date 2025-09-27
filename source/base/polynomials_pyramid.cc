// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/base/polynomials_barycentric.h>
#include <deal.II/base/polynomials_pyramid.h>

DEAL_II_NAMESPACE_OPEN

namespace
{
  unsigned int
  compute_n_polynomials_pyramid(const unsigned int dim,
                                const unsigned int degree)
  {
    if (dim == 3)
      {
        if (degree == 1)
          return 5;
      }

    DEAL_II_NOT_IMPLEMENTED();

    return 0;
  }
} // namespace



template <int dim>
ScalarLagrangePolynomialPyramid<dim>::ScalarLagrangePolynomialPyramid(
  const unsigned int degree)
  : ScalarPolynomialsBase<dim>(degree,
                               compute_n_polynomials_pyramid(dim, degree))
{}


template <int dim>
double
ScalarLagrangePolynomialPyramid<dim>::compute_value(const unsigned int i,
                                                    const Point<dim>  &p) const
{
  AssertDimension(dim, 3);
  AssertIndexRange(this->degree(), 2);

  const double Q14 = 0.25;
  double       ration;

  const double r = p[0];
  const double s = p[1];
  const double t = p[2];

  if (std::fabs(t - 1.0) > 1.0e-14)
    {
      ration = (r * s * t) / (1.0 - t);
    }
  else
    {
      ration = 0.0;
    }

  if (i == 0)
    return Q14 * ((1.0 - r) * (1.0 - s) - t + ration);
  if (i == 1)
    return Q14 * ((1.0 + r) * (1.0 - s) - t - ration);
  if (i == 2)
    return Q14 * ((1.0 - r) * (1.0 + s) - t - ration);
  if (i == 3)
    return Q14 * ((1.0 + r) * (1.0 + s) - t + ration);
  else
    return t;
}



template <int dim>
Tensor<1, dim>
ScalarLagrangePolynomialPyramid<dim>::compute_grad(const unsigned int i,
                                                   const Point<dim>  &p) const
{
  AssertDimension(dim, 3);
  AssertIndexRange(this->degree(), 4);

  Tensor<1, dim> grad;

  if (this->degree() == 1)
    {
      const double Q14 = 0.25;

      const double r = p[0];
      const double s = p[1];
      const double t = p[2];

      double rationdr;
      double rationds;
      double rationdt;

      if (std::fabs(t - 1.0) > 1.0e-14)
        {
          rationdr = s * t / (1.0 - t);
          rationds = r * t / (1.0 - t);
          rationdt = r * s / ((1.0 - t) * (1.0 - t));
        }
      else
        {
          rationdr = 1.0;
          rationds = 1.0;
          rationdt = 1.0;
        }


      if (i == 0)
        {
          grad[0] = Q14 * (-1.0 * (1.0 - s) + rationdr);
          grad[1] = Q14 * (-1.0 * (1.0 - r) + rationds);
          grad[2] = Q14 * (rationdt - 1.0);
        }
      else if (i == 1)
        {
          grad[0] = Q14 * (1.0 * (1.0 - s) - rationdr);
          grad[1] = Q14 * (-1.0 * (1.0 + r) - rationds);
          grad[2] = Q14 * (-1.0 * rationdt - 1.0);
        }
      else if (i == 2)
        {
          grad[0] = Q14 * (-1.0 * (1.0 + s) - rationdr);
          grad[1] = Q14 * (1.0 * (1.0 - r) - rationds);
          grad[2] = Q14 * (-1.0 * rationdt - 1.0);
        }
      else if (i == 3)
        {
          grad[0] = Q14 * (1.0 * (1.0 + s) + rationdr);
          grad[1] = Q14 * (1.0 * (1.0 + r) + rationds);
          grad[2] = Q14 * (rationdt - 1.0);
        }
      else if (i == 4)
        {
          grad[0] = 0.0;
          grad[1] = 0.0;
          grad[2] = 1.0;
        }
      else
        {
          DEAL_II_NOT_IMPLEMENTED();
        }
    }

  return grad;
}



template <int dim>
Tensor<2, dim>
ScalarLagrangePolynomialPyramid<dim>::compute_grad_grad(
  const unsigned int i,
  const Point<dim>  &p) const
{
  (void)i;
  (void)p;

  DEAL_II_NOT_IMPLEMENTED();
  return Tensor<2, dim>();
}



template <int dim>
void
ScalarLagrangePolynomialPyramid<dim>::evaluate(
  const Point<dim>            &unit_point,
  std::vector<double>         &values,
  std::vector<Tensor<1, dim>> &grads,
  std::vector<Tensor<2, dim>> &grad_grads,
  std::vector<Tensor<3, dim>> &third_derivatives,
  std::vector<Tensor<4, dim>> &fourth_derivatives) const
{
  (void)grads;
  (void)grad_grads;
  (void)third_derivatives;
  (void)fourth_derivatives;

  if (values.size() == this->n())
    for (unsigned int i = 0; i < this->n(); ++i)
      values[i] = compute_value(i, unit_point);

  if (grads.size() == this->n())
    for (unsigned int i = 0; i < this->n(); ++i)
      grads[i] = compute_grad(i, unit_point);
}



template <int dim>
Tensor<1, dim>
ScalarLagrangePolynomialPyramid<dim>::compute_1st_derivative(
  const unsigned int i,
  const Point<dim>  &p) const
{
  return compute_grad(i, p);
}



template <int dim>
Tensor<2, dim>
ScalarLagrangePolynomialPyramid<dim>::compute_2nd_derivative(
  const unsigned int i,
  const Point<dim>  &p) const
{
  (void)i;
  (void)p;

  DEAL_II_NOT_IMPLEMENTED();

  return {};
}



template <int dim>
Tensor<3, dim>
ScalarLagrangePolynomialPyramid<dim>::compute_3rd_derivative(
  const unsigned int i,
  const Point<dim>  &p) const
{
  (void)i;
  (void)p;

  DEAL_II_NOT_IMPLEMENTED();

  return {};
}



template <int dim>
Tensor<4, dim>
ScalarLagrangePolynomialPyramid<dim>::compute_4th_derivative(
  const unsigned int i,
  const Point<dim>  &p) const
{
  (void)i;
  (void)p;

  DEAL_II_NOT_IMPLEMENTED();

  return {};
}



template <int dim>
std::string
ScalarLagrangePolynomialPyramid<dim>::name() const
{
  return "ScalarLagrangePolynomialPyramid";
}



template <int dim>
std::unique_ptr<ScalarPolynomialsBase<dim>>
ScalarLagrangePolynomialPyramid<dim>::clone() const
{
  return std::make_unique<ScalarLagrangePolynomialPyramid<dim>>(*this);
}



template class ScalarLagrangePolynomialPyramid<1>;
template class ScalarLagrangePolynomialPyramid<2>;
template class ScalarLagrangePolynomialPyramid<3>;

DEAL_II_NAMESPACE_CLOSE
