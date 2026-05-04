// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2020 - 2025 by the deal.II authors
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
#include <deal.II/base/polynomials_barycentric.h>
#include <deal.II/base/polynomials_pyramid.h>
#include <deal.II/base/scalar_polynomials_vandermonde_base.h>
#include <deal.II/base/table.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/utilities.h>

#include <deal.II/grid/reference_cell.h>

#include <Kokkos_Macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <memory>
#include <string>
#include <vector>

DEAL_II_NAMESPACE_OPEN


namespace
{
  /**
   * Return the vertices of a pyramid in a way that can also be compiled
   * for dim==1 and dim==2.
   */
  template <int dim>
  std::vector<Point<dim>>
  get_pyramid_vertices()
  {
    if constexpr (dim == 3)
      return {ReferenceCells::Pyramid.vertex(0),
              ReferenceCells::Pyramid.vertex(1),
              ReferenceCells::Pyramid.vertex(2),
              ReferenceCells::Pyramid.vertex(3),
              ReferenceCells::Pyramid.vertex(4)};
    else
      {
        DEAL_II_ASSERT_UNREACHABLE();
        return {};
      }
  }
} // namespace



template <int dim>
ScalarLagrangePolynomialPyramid<dim>::ScalarLagrangePolynomialPyramid(
  const unsigned int degree)
  : ScalarLagrangePolynomialPyramid<dim>(1, 5, get_pyramid_vertices<dim>())
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
  : ScalarPolynomialsVandermondeBase<dim>(degree, n_dofs)
{
  AssertThrow(dim == 3,
              ExcNotImplemented("Pyramid elements only make sense in 3D."));

  this->reinit(support_points);
}



template <int dim>
double
ScalarLagrangePolynomialPyramid<
  dim>::evaluate_orthogonal_basis_function_by_degree(const unsigned int i,
                                                     const unsigned int j,
                                                     const unsigned int k,
                                                     const Point<dim>  &p) const
{
  AssertIndexRange(i, this->degree() + 1);
  AssertIndexRange(j, this->degree() + 1);
  AssertIndexRange(k, this->degree() + 1 - std::max(i, j));

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
ScalarLagrangePolynomialPyramid<dim>::evaluate_orthogonal_basis_function(
  const unsigned int i,
  const Point<dim>  &p) const
{
  AssertIndexRange(i, this->n());

  // find corresponding entry to i
  // it holds 0 <= j, k <= degree,
  // 0 <= l <= degree - max(j, k)
  for (unsigned int j = 0, counter = 0; j <= this->degree(); ++j)
    for (unsigned int k = 0; k <= this->degree(); ++k)
      for (unsigned int l = 0; l <= this->degree() - std::max(j, k);
           ++l, ++counter)
        if (counter == i)
          return evaluate_orthogonal_basis_function_by_degree(j, k, l, p);

  DEAL_II_ASSERT_UNREACHABLE();
  return 0;
}



template <int dim>
Tensor<1, dim>
ScalarLagrangePolynomialPyramid<dim>::
  evaluate_orthogonal_basis_derivative_by_degree(const unsigned int i,
                                                 const unsigned int j,
                                                 const unsigned int k,
                                                 const Point<dim>  &p) const
{
  AssertIndexRange(i, this->degree() + 1);
  AssertIndexRange(j, this->degree() + 1);
  AssertIndexRange(k, this->degree() + 1 - std::max(i, j));

  Tensor<1, dim> grad;

  const double x = p[0];
  const double y = p[1];
  const double z = p[2];

  // handle the special cases where 1/(1-z) cancels and the tip
  if (i == 0 && j == 0)
    {
      // assume (1-z)^0 = 1 even when z = 1
      grad[0] = 0.;
      if constexpr (dim > 1)
        grad[1] = 0.;
      if constexpr (dim > 2)
        grad[2] =
          Polynomials::jacobi_polynomial_derivative<double>(k, 2, 0, z, true);

      return grad;
    }
  else if (i == 0 && j == 1)
    {
      grad[0] = 0.;
      if constexpr (dim > 1)
        grad[1] =
          Polynomials::jacobi_polynomial_value<double>(k, 4, 0, z, true);
      if constexpr (dim > 2)
        grad[2] =
          y *
          Polynomials::jacobi_polynomial_derivative<double>(k, 4, 0, z, true);

      return grad;
    }
  else if (i == 1 && j == 0)
    {
      grad[0] = Polynomials::jacobi_polynomial_value<double>(k, 4, 0, z, true);
      if constexpr (dim > 1)
        grad[1] = 0.0;
      if constexpr (dim > 2)
        grad[2] =
          x *
          Polynomials::jacobi_polynomial_derivative<double>(k, 4, 0, z, true);

      return grad;
    }
  else if (std::abs(p[2] - 1.0) < 1e-14)
    {
      grad = 0.;
      if (i == 1 && j == 1)
        // assume x/(1-z)-> 1 and y/(1-z)-> 1 at the tip
        for (unsigned int d = 0; d < dim; ++d)
          grad[d] =
            Polynomials::jacobi_polynomial_value<double>(k, 4, 0, z, true);

      return grad;
    }

  const unsigned int max_ij = std::max(i, j);
  const double       ratio  = 1.0 / (1.0 - z);

  grad[0] =
    Polynomials::jacobi_polynomial_derivative<double>(
      i, 0, 0, x * ratio, false) *
    ratio *
    Polynomials::jacobi_polynomial_value<double>(j, 0, 0, y * ratio, false) *
    std::pow((1.0 - z), max_ij) *
    Polynomials::jacobi_polynomial_value<double>(k, 2 * max_ij + 2, 0, z, true);
  if constexpr (dim > 1)
    grad[1] =
      Polynomials::jacobi_polynomial_value<double>(i, 0, 0, x * ratio, false) *
      Polynomials::jacobi_polynomial_derivative<double>(
        j, 0, 0, y * ratio, false) *
      ratio * std::pow((1.0 - z), max_ij) *
      Polynomials::jacobi_polynomial_value<double>(
        k, 2 * max_ij + 2, 0, z, true);
  if constexpr (dim > 2)
    grad[2] =
      Polynomials::jacobi_polynomial_derivative<double>(
        i, 0, 0, x * ratio, false) *
        x * ratio * ratio *
        Polynomials::jacobi_polynomial_value<double>(
          j, 0, 0, y * ratio, false) *
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
          k, 2 * max_ij + 2, 0, 2.0 * z - 1.0, false) *
        2.0;

  for (unsigned int d = 0; d < dim; ++d)
    if (std::fabs(grad[d]) < 1e-14)
      grad[d] = 0.0;

  return grad;
}



template <int dim>
Tensor<1, dim>
ScalarLagrangePolynomialPyramid<dim>::evaluate_orthogonal_basis_derivative(
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
          return evaluate_orthogonal_basis_derivative_by_degree(j, k, l, p);

  DEAL_II_ASSERT_UNREACHABLE();
  return Tensor<1, dim>();
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
