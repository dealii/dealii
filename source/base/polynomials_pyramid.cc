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


#include <deal.II/base/point.h>
#include <deal.II/base/polynomial.h>
#include <deal.II/base/polynomials_barycentric.h>
#include <deal.II/base/polynomials_pyramid.h>

#include <deal.II/grid/reference_cell.h>

DEAL_II_NAMESPACE_OPEN



template <int dim>
ScalarLagrangePolynomialPyramid<dim>::ScalarLagrangePolynomialPyramid(
  const unsigned int degree)
  : ScalarLagrangePolynomialPyramid<dim>(
      1,
      5,
      {ReferenceCells::Pyramid.vertex<dim>(0),
       ReferenceCells::Pyramid.vertex<dim>(1),
       ReferenceCells::Pyramid.vertex<dim>(2),
       ReferenceCells::Pyramid.vertex<dim>(3),
       ReferenceCells::Pyramid.vertex<dim>(4)})
{
  AssertThrow(degree == 1,
              ExcNotImplemented(
                "This constructor only works for linear elements."));
}



template <int dim>
ScalarLagrangePolynomialPyramid<dim>::ScalarLagrangePolynomialPyramid(
  const unsigned int             degree,
  const unsigned int             n_dofs,
  const std::vector<Point<dim>> &support_points)
  : ScalarPolynomialsBase<dim>(degree, n_dofs)
{
  (void)n_dofs;

  AssertThrow(dim == 3,
              ExcNotImplemented("Pyramid elements only make sense in 3D."));

  AssertThrow(
    support_points.size() == n_dofs,
    ExcNotImplemented(
      "The number of DoF must be equal to the number of support points."));

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

  vandermonde_matrix_inverse = VDM;
}



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

  const unsigned int max_ij = std::max(i, j);

  // handle the special cases
  // assume (1-z)^0 = 1 even when z = 1
  if (max_ij == 0)
    return Polynomials::jacobi_polynomial_value<double>(k, 2, 0, z, true);

  // at the tip the basis looks like (x/(1-z))^i * (y/(1-z))^j * (1-z)^max(i,j)
  // when going to the tip of the cone the relation is roughly x = 1 - z, so
  // (x/(1-z)) = const and the (1-z)^max(i,j) term then pushes to 0
  if (std::fabs(z - 1.0) < 1e-14)
    return 0.0;

  const double ratio = 1.0 / (1.0 - z);

  const double phi =
    Polynomials::jacobi_polynomial_value<double>(i, 0, 0, x * ratio, false) *
    Polynomials::jacobi_polynomial_value<double>(j, 0, 0, y * ratio, false) *
    std::pow((1.0 - z), max_ij) *
    Polynomials::jacobi_polynomial_value<double>(k, 2 * max_ij + 2, 0, z, true);

  if (std::fabs(phi) < 1e-14)
    return 0.0;

  return phi;
}



template <int dim>
double
ScalarLagrangePolynomialPyramid<dim>::compute_jacobi_basis(
  const unsigned int i,
  const Point<dim>  &p) const
{
  AssertIndexRange(i, this->n());

  // find corresponding entry to i
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
ScalarLagrangePolynomialPyramid<dim>::compute_polynomial_space_derivative(
  const unsigned int i,
  const unsigned int j,
  const unsigned int k,
  const Point<dim>  &p) const
{
  static_assert(dim == 3,
    "Pyramid elements only make sense for dim=3.");

  AssertThrow(std::abs(p[2] - 1.0) > 1e-14,
              ExcMessage("The derivative at the tip is not defined."));
  // e.g. for i,j,k=1,1,0 the gradient is [y/(1-z),x/(1-z),xy/(1-z)**2] which is
  // not defined at [0,0,1]

  const double x = p[0];
  const double y = p[1];
  const double z = p[2];

  Tensor<1, dim> grad;

  const unsigned int max_ij = std::max(i, j);
  const double       ratio  = 1.0 / (1.0 - z);

  grad[0] =
    Polynomials::jacobi_polynomial_derivative<double>(
      i, 0, 0, x * ratio, false) *
    ratio *
    Polynomials::jacobi_polynomial_value<double>(j, 0, 0, y * ratio, false) *
    std::pow((1.0 - z), max_ij) *
    Polynomials::jacobi_polynomial_value<double>(k, 2 * max_ij + 2, 0, z, true);
  grad[1] =
    Polynomials::jacobi_polynomial_value<double>(i, 0, 0, x * ratio, false) *
    Polynomials::jacobi_polynomial_derivative<double>(
      j, 0, 0, y * ratio, false) *
    ratio * std::pow((1.0 - z), max_ij) *
    Polynomials::jacobi_polynomial_value<double>(k, 2 * max_ij + 2, 0, z, true);
  grad[2] =
    Polynomials::jacobi_polynomial_derivative<double>(
      i, 0, 0, x * ratio, false) *
      x * ratio * ratio *
      Polynomials::jacobi_polynomial_value<double>(j, 0, 0, y * ratio, false) *
      std::pow((1.0 - z), max_ij) *
      Polynomials::jacobi_polynomial_value<double>(
        k, 2 * max_ij + 2, 0, z, true) +
    Polynomials::jacobi_polynomial_value<double>(i, 0, 0, x * ratio, false) *
      Polynomials::jacobi_polynomial_derivative<double>(
        j, 0, 0, y * ratio, false) *
      y * ratio * ratio * std::pow((1.0 - z), max_ij) *
      Polynomials::jacobi_polynomial_value<double>(
        k, 2 * max_ij + 2, 0, z, true) +
    Polynomials::jacobi_polynomial_value<double>(i, 0, 0, x * ratio, false) *
      Polynomials::jacobi_polynomial_value<double>(j, 0, 0, y * ratio, false) *
      (-1.0) * max_ij * std::pow((1.0 - z), max_ij - 1) *
      Polynomials::jacobi_polynomial_value<double>(
        k, 2 * max_ij + 2, 0, z, true) +
    Polynomials::jacobi_polynomial_value<double>(i, 0, 0, x * ratio, false) *
      Polynomials::jacobi_polynomial_value<double>(j, 0, 0, y * ratio, false) *
      std::pow((1.0 - z), max_ij) *
      Polynomials::jacobi_polynomial_derivative<double>(
        k, 2 * max_ij + 2, 0, 2.0 * z - 1.0, false) *
      2.0;

  for (unsigned int d = 0; d < dim; ++d)
    if (std::fabs(grad[d]) < 1e-14)
      grad[d] = 0.0;

  return grad;
}



template <int dim>
Tensor<1, dim>
ScalarLagrangePolynomialPyramid<dim>::compute_jacobi_derivative(
  const unsigned int i,
  const Point<dim>  &p) const
{
  AssertIndexRange(i, this->n());

  // find corresponding entrance to i
  for (unsigned int j = 0, counter = 0; j <= this->degree(); ++j)
    for (unsigned int k = 0; k <= this->degree(); ++k)
      for (unsigned int l = 0; l <= this->degree() - std::max(j, k);
           ++l, ++counter)
        if (counter == i)
          return compute_polynomial_space_derivative(j, k, l, p);

  DEAL_II_ASSERT_UNREACHABLE();
  return Tensor<1, dim>();
}



template <int dim>
double
ScalarLagrangePolynomialPyramid<dim>::compute_value(const unsigned int i,
                                                    const Point<dim>  &p) const
{
  AssertDimension(dim, 3);
  AssertIndexRange(i, vandermonde_matrix_inverse.m());

  double result = 0;
  for (unsigned int j = 0; j < vandermonde_matrix_inverse.n(); ++j)
    result +=
      vandermonde_matrix_inverse[i][j] * this->compute_jacobi_basis(j, p);

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
  AssertIndexRange(i, vandermonde_matrix_inverse.m());

  Tensor<1, dim> grad;

  for (unsigned int j = 0; j < vandermonde_matrix_inverse.n(); ++j)
    grad +=
      vandermonde_matrix_inverse[i][j] * this->compute_jacobi_derivative(j, p);

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
