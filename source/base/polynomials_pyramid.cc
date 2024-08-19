// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/base/point.h>
#include <deal.II/base/polynomial.h>
#include <deal.II/base/polynomials_barycentric.h>
#include <deal.II/base/polynomials_pyramid.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/reference_cell.h>


DEAL_II_NAMESPACE_OPEN


template <int dim>
double
ScalarLagrangePolynomialPyramid<dim>::compute_polynomial_space(
  const unsigned int i,
  const unsigned int j,
  const unsigned int k,
  const Point<dim>  &p) const
{
  const double x = p[0];
  const double y = p[1];
  const double z = p[2];

  double ratio;
  if (std::fabs(z - 1.0) < 1e-12)
    ratio = 0.0;
  else
    ratio = 1.0 / (1.0 - z);

  double             phi    = 0.0;
  const unsigned int max_ij = std::max(i, j);

  phi =
    Polynomials::jacobi_polynomial_value<double>(i, 0, 0, x * ratio, false) *
    Polynomials::jacobi_polynomial_value<double>(j, 0, 0, y * ratio, false) *
    std::pow((1.0 - z), max_ij) *
    Polynomials::jacobi_polynomial_value<double>(k, 2 * max_ij + 2, 0, z, true);

  if (std::fabs(phi) < 1e-12)
    return 0.0;

  return phi;
}



template <int dim>
double
ScalarLagrangePolynomialPyramid<dim>::compute_jacobi_basis(
  const unsigned int i,
  const Point<dim>  &p) const
{
  AssertIndexRange(i,
                   (this->degree() + 1) * (this->degree() + 2) *
                     (2 * this->degree() + 3) / 6);

  // find corresponding entrance to i
  for (unsigned int j = 0, counter = 0; j <= this->degree(); ++j)
    for (unsigned int k = 0; k <= this->degree(); ++k)
      for (unsigned int l = 0; l <= this->degree() - std::max(j, k);
           ++l, ++counter)
        if (counter == i)
          return compute_polynomial_space(j, k, l, p);

  DEAL_II_ASSERT_UNREACHABLE();
  return 0;
}



template <int dim>
Tensor<1, dim>
ScalarLagrangePolynomialPyramid<dim>::compute_polynomial_space_deriv(
  const unsigned int i,
  const unsigned int j,
  const unsigned int k,
  const Point<dim>  &p) const
{
  const double x = p[0];
  const double y = p[1];
  const double z = p[2];

  Tensor<1, dim> grad;
  double         ratio;
  if (std::fabs(z - 1.0) < 1e-12)
    ratio = 0.0;
  else
    ratio = 1.0 / (1.0 - z);

  const unsigned int max_ij = std::max(i, j);

  grad[0] =
    Polynomials::jacobi_polynomial_derivative<double>(
      i, 0, 0, x * ratio, false) *
    ratio *
    Polynomials::jacobi_polynomial_value<double>(j, 0, 0, y * ratio, false) *
    std::pow((1.0 - z), max_ij) *
    Polynomials::jacobi_polynomial_value<double>(k, 2 * max_ij + 2, 0, z, true);
  grad[1] =
    Polynomials::jacobi_polynomial_derivative<double>(
      j, 0, 0, y * ratio, false) *
    ratio *
    Polynomials::jacobi_polynomial_value<double>(i, 0, 0, x * ratio, false) *
    std::pow((1.0 - z), max_ij) *
    Polynomials::jacobi_polynomial_value<double>(k, 2 * max_ij + 2, 0, z, true);
  if (max_ij == 0)
    grad[2] =
      -x * std::pow(ratio, 2) *
        Polynomials::jacobi_polynomial_derivative<double>(
          i, 0, 0, x * ratio, false) *
        Polynomials::jacobi_polynomial_value<double>(
          j, 0, 0, y * ratio, false) *
        Polynomials::jacobi_polynomial_value<double>(k, 2, 0, z, true) -
      y * std::pow(ratio, 2) *
        Polynomials::jacobi_polynomial_derivative<double>(
          j, 0, 0, y * ratio, false) *
        Polynomials::jacobi_polynomial_value<double>(
          i, 0, 0, x * ratio, false) *
        Polynomials::jacobi_polynomial_value<double>(k, 2, 0, z, true) +
      Polynomials::jacobi_polynomial_value<double>(i, 0, 0, x * ratio, false) *
        Polynomials::jacobi_polynomial_value<double>(
          j, 0, 0, y * ratio, false) *
        Polynomials::jacobi_polynomial_derivative<double>(
          k, 2 * max_ij + 2, 0, z, true);
  else
    grad[2] =
      -x * std::pow(ratio, 2) *
        Polynomials::jacobi_polynomial_derivative<double>(
          i, 0, 0, x * ratio, false) *
        Polynomials::jacobi_polynomial_value<double>(
          j, 0, 0, y * ratio, false) *
        std::pow((1.0 - z), max_ij) *
        Polynomials::jacobi_polynomial_value<double>(k, 2, 0, z, true) -
      y * std::pow(ratio, 2) *
        Polynomials::jacobi_polynomial_derivative<double>(
          j, 0, 0, y * ratio, false) *
        Polynomials::jacobi_polynomial_value<double>(
          i, 0, 0, x * ratio, false) *
        std::pow((1.0 - z), max_ij) *
        Polynomials::jacobi_polynomial_value<double>(k, 2, 0, z, true) +
      Polynomials::jacobi_polynomial_value<double>(i, 0, 0, x * ratio, false) *
        Polynomials::jacobi_polynomial_value<double>(
          j, 0, 0, y * ratio, false) *
        (-1.0) * max_ij * std::pow((1.0 - z), max_ij - 1) *
        Polynomials::jacobi_polynomial_value<double>(
          k, 2 * max_ij + 2, 0, z, true) +
      Polynomials::jacobi_polynomial_value<double>(i, 0, 0, x * ratio, false) *
        Polynomials::jacobi_polynomial_value<double>(
          j, 0, 0, y * ratio, false) *
        std::pow((1.0 - z), max_ij) *
        Polynomials::jacobi_polynomial_derivative<double>(
          k, 2 * max_ij + 2, 0, z, true);


  for (unsigned int d = 0; d < dim; ++d)
    if (std::fabs(grad[d]) < 1e-12)
      grad[d] = 0.0;

  return grad;
}



template <int dim>
Tensor<1, dim>
ScalarLagrangePolynomialPyramid<dim>::compute_jacobi_deriv(
  const unsigned int i,
  const Point<dim>  &p) const
{
  AssertIndexRange(i,
                   (this->degree() + 1) * (this->degree() + 2) *
                     (2 * this->degree() + 3) / 6);

  // find corresponding entrance to i
  for (unsigned int j = 0, counter = 0; j <= this->degree(); ++j)
    for (unsigned int k = 0; k <= this->degree(); ++k)
      for (unsigned int l = 0; l <= this->degree() - std::max(j, k);
           ++l, ++counter)
        if (counter == i)
          return compute_polynomial_space_deriv(j, k, l, p);

  DEAL_II_ASSERT_UNREACHABLE();
  return Tensor<1, dim>();
}



template <int dim>
ScalarLagrangePolynomialPyramid<dim>::ScalarLagrangePolynomialPyramid(
  const unsigned int            degree,
  const unsigned int            n_dofs,
  const std::vector<Point<dim>> support_points)
  : ScalarPolynomialsBase<dim>(degree, n_dofs)
{
  // fill VDM Matrix
  FullMatrix<double> VDM(support_points.size());

  for (unsigned i = 0; i < VDM.m(); ++i)
    for (unsigned j = 0; j < VDM.n(); ++j)
      VDM[i][j] = this->compute_jacobi_basis(i, support_points[j]);

  // get inverse
  VDM.gauss_jordan();

  for (unsigned i = 0; i < VDM.m(); ++i)
    for (unsigned j = 0; j < VDM.n(); ++j)
      if (std::fabs(VDM[i][j]) < 1e-14)
        VDM[i][j] = 0.0;

  VDM_inv = VDM;
}



template <int dim>
double
ScalarLagrangePolynomialPyramid<dim>::compute_value(const unsigned int i,
                                                    const Point<dim>  &p) const
{
  AssertDimension(dim, 3);
  AssertIndexRange(i, VDM_inv.m());

  double result = 0;
  for (unsigned int j = 0; j < VDM_inv.n(); ++j)
    result += VDM_inv[i][j] * this->compute_jacobi_basis(j, p);

  if (std::fabs(result) < 1e-14)
    result = 0.0;

  return result;
}



template <int dim>
Tensor<1, dim>
ScalarLagrangePolynomialPyramid<dim>::compute_grad(const unsigned int i,
                                                   const Point<dim>  &p) const
{
  AssertDimension(dim, 3);
  AssertIndexRange(i, VDM_inv.m());

  Tensor<1, dim> grad;

  for (unsigned int j = 0; j < VDM_inv.n(); ++j)
    grad += VDM_inv[i][j] * this->compute_jacobi_deriv(j, p);

  for (unsigned int d = 0; d < dim; ++d)
    if (std::fabs(grad[d]) < 1e-14)
      grad[d] = 0.0;

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
