// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2020 - 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


#include <deal.II/base/config.h>

#include <deal.II/base/exception_macros.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/numbers.h>
#include <deal.II/base/point.h>
#include <deal.II/base/polynomial.h>
#include <deal.II/base/polynomials_wedge.h>
#include <deal.II/base/scalar_polynomials_base.h>
#include <deal.II/base/scalar_polynomials_vandermonde_base.h>
#include <deal.II/base/table.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/utilities.h>

#include <Kokkos_Macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <memory>
#include <string>
#include <vector>

DEAL_II_NAMESPACE_OPEN



template <int dim>
ScalarLagrangePolynomialWedge<dim>::ScalarLagrangePolynomialWedge(
  const unsigned int             degree,
  const std::vector<Point<dim>> &support_points)
  : ScalarPolynomialsVandermondeBase<dim>(degree, support_points.size())
{
  AssertDimension(dim, 3);

  const unsigned int n_dofs = (degree + 1) * (degree + 1) * (degree + 2) / 2;
  AssertDimension(n_dofs, support_points.size());

  this->reinit(support_points);
}



template <int dim>
ScalarLagrangePolynomialWedge<dim>::ScalarLagrangePolynomialWedge(
  const unsigned int degree)
  : ScalarLagrangePolynomialWedge<dim>(degree,
                                       internal::get_wedge_support_points<dim>(
                                         degree))
{
  AssertThrow(
    degree == 1 || degree == 2,
    ExcNotImplemented(
      "This constructor only works for linear and quadratic elements."));
}



template <int dim>
double
ScalarLagrangePolynomialWedge<
  dim>::evaluate_orthogonal_basis_function_by_degree(const unsigned int i,
                                                     const unsigned int j,
                                                     const unsigned int k,
                                                     const Point<dim>  &p) const
{
  AssertIndexRange(i + j, this->degree() + 1);
  AssertIndexRange(k, this->degree() + 1);

  const double x = p[0];
  const double y = p[1];
  const double z = p[2];

  // the basis function looks like
  // P_i^{0,0}(2x/(1-y)-1) * (1-y)^i * P_j^{2*i+1,0}(2*y-1) * P_k^{0,0}(2z-1)
  // separate it into
  // P_i^{0,0}(2x/(1-y)-1) * (1-y)^i
  // P_j^{2*i+1,0}(2*y-1)
  // P_k^{0,0}(2z-1)
  // the first is a homogenized Jacobi polynomial
  // define s = 1 - y so the first term can be written as
  // Q_i^{0,0}(x,s) = P_i^{0,0}(2 * x/s - 1) * s^i
  const double s = 1 - y;

  const double Qi =
    Polynomials::jacobi_polynomial_homogenized_value<double>(i, 0, 0, x, s);

  const double Pj =
    Polynomials::jacobi_polynomial_value<double>(j, 2 * i + 1, 0, y, true);

  const double Pk =
    Polynomials::jacobi_polynomial_value<double>(k, 0, 0, z, true);

  const double phi = Qi * Pj * Pk;

  if (std::fabs(phi) < 1e-14)
    return 0.0;

  return phi;
}



template <int dim>
double
ScalarLagrangePolynomialWedge<dim>::evaluate_orthogonal_basis_function(
  const unsigned int i,
  const Point<dim>  &p) const
{
  AssertIndexRange(i, this->n());

  // find corresponding entry to i
  // it holds 0 <= j + k <= degree
  // 0 <= l <= degree
  for (unsigned int j = 0, counter = 0; j < this->degree() + 1; ++j)
    for (unsigned int k = 0; k < this->degree() + 1 - j; ++k)
      for (unsigned int l = 0; l < this->degree() + 1; ++l, ++counter)
        if (counter == i)
          return evaluate_orthogonal_basis_function_by_degree(j, k, l, p);

  DEAL_II_ASSERT_UNREACHABLE();
  return 0;
}



template <int dim>
Tensor<1, dim>
ScalarLagrangePolynomialWedge<dim>::
  evaluate_orthogonal_basis_derivative_by_degree(const unsigned int i,
                                                 const unsigned int j,
                                                 const unsigned int k,
                                                 const Point<dim>  &p) const
{
  AssertIndexRange(i + j, this->degree() + 1);
  AssertIndexRange(k, this->degree() + 1);

  Tensor<1, dim> grad;

  const double x = p[0];
  const double y = p[1];
  const double z = p[2];

  // the basis function looks like
  // P_i^{0,0}(2x/(1-y)-1) * (1-y)^i * P_j^{2*i+1,0}(2*y-1) * P_k^{0,0}(2z-1)
  // separate it into
  // P_i^{0,0}(2x/(1-y)-1) * (1-y)^i
  // and
  // P_j^{2*i+1,0}(2*y-1)
  // the first is a homogenized Jacobi polynomial
  // define s = 1 - y so the first term can be written as
  // Q_i^{0,0}(x,s) = P_i^{0,0}(2 * x/s - 1) * s^i

  // to get the derivatives just use the product rule with all terms
  const double s     = 1 - y;
  const double ds_dy = -1.0;

  const double Qi =
    Polynomials::jacobi_polynomial_homogenized_value<double>(i, 0, 0, x, s);
  const double Pj =
    Polynomials::jacobi_polynomial_value<double>(j, 2 * i + 1, 0, y, true);
  const double Pk =
    Polynomials::jacobi_polynomial_value<double>(k, 0, 0, z, true);

  const auto dQi_dx =
    Polynomials::jacobi_polynomial_homogenized_derivative<double>(
      1, 0, i, 0, 0, x, s);

  const auto dQi_ds =
    Polynomials::jacobi_polynomial_homogenized_derivative<double>(
      0, 1, i, 0, 0, x, s);

  const auto dPj_dy =
    Polynomials::jacobi_polynomial_derivative<double>(j, 2 * i + 1, 0, y, true);

  const double dPk_dz =
    Polynomials::jacobi_polynomial_derivative<double>(k, 0, 0, z, true);

  grad[0] = dQi_dx * Pj * Pk;
  if constexpr (dim > 1)
    grad[1] = dQi_ds * ds_dy * Pj * Pk + Qi * dPj_dy * Pk;
  if constexpr (dim > 2)
    grad[2] = Qi * Pj * dPk_dz;

  for (unsigned int d = 0; d < dim; ++d)
    if (std::fabs(grad[d]) < 1e-14)
      grad[d] = 0.0;

  return grad;
}



template <int dim>
Tensor<1, dim>
ScalarLagrangePolynomialWedge<dim>::evaluate_orthogonal_basis_derivative(
  const unsigned int i,
  const Point<dim>  &p) const
{
  AssertIndexRange(i, this->n());

  // find corresponding entry to i
  // it holds 0 <= j + k <= degree
  // 0 <= l <= degree
  for (unsigned int j = 0, counter = 0; j < this->degree() + 1; ++j)
    for (unsigned int k = 0; k < this->degree() + 1 - j; ++k)
      for (unsigned int l = 0; l < this->degree() + 1; ++l, ++counter)
        if (counter == i)
          if (counter == i)
            return evaluate_orthogonal_basis_derivative_by_degree(j, k, l, p);

  DEAL_II_ASSERT_UNREACHABLE();
  return Tensor<1, dim>();
}



template <int dim>
std::string
ScalarLagrangePolynomialWedge<dim>::name() const
{
  return "ScalarLagrangePolynomialWedge";
}



template <int dim>
std::unique_ptr<ScalarPolynomialsBase<dim>>
ScalarLagrangePolynomialWedge<dim>::clone() const
{
  return std::make_unique<ScalarLagrangePolynomialWedge<dim>>(*this);
}



template class ScalarLagrangePolynomialWedge<1>;
template class ScalarLagrangePolynomialWedge<2>;
template class ScalarLagrangePolynomialWedge<3>;

DEAL_II_NAMESPACE_CLOSE
