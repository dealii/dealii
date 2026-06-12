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
#include <deal.II/base/polynomials_barycentric.h>
#include <deal.II/base/polynomials_simplex.h>
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


template <int dim>
ScalarLagrangePolynomialSimplex<dim>::ScalarLagrangePolynomialSimplex(
  const unsigned int             degree,
  const std::vector<Point<dim>> &support_points)
  : ScalarPolynomialsVandermondeBase<dim>(degree, support_points.size())
{
  if constexpr (dim == 1)
    AssertDimension(degree + 1, support_points.size());
  else if constexpr (dim == 2)
    {
      const unsigned int n_dofs = (degree + 1) * (degree + 2) / 2;
      AssertDimension(n_dofs, support_points.size());
    }
  else if constexpr (dim == 3)
    {
      const unsigned int n_dofs =
        (degree + 1) * (degree + 2) * (degree + 3) / 6;
      AssertDimension(n_dofs, support_points.size());
    }
  else
    DEAL_II_NOT_IMPLEMENTED();

  this->reinit(support_points);
}



template <int dim>
double
ScalarLagrangePolynomialSimplex<
  dim>::evaluate_orthogonal_basis_function_by_degree(const unsigned int i,
                                                     const unsigned int j,
                                                     const unsigned int k,
                                                     const Point<dim>  &p) const
{
  AssertIndexRange(i + j + k, this->degree() + 1);

  const double x = p[0];

  if constexpr (dim == 1)
    {
      const double phi =
        Polynomials::jacobi_polynomial_value<double>(i, 0, 0, x, true);

      if (std::fabs(phi) < 1e-14)
        return 0.0;

      return phi;
    }
  else if constexpr (dim == 2)
    {
      const double y      = p[1];
      const double factor = std::abs(1.0 - y) < 1e-14 ? 1.0 : 1.0 / (1.0 - y);

      const double x_contribution =
        i == 0 ? 1.0 :
                 Polynomials::jacobi_polynomial_value<double>(
                   i, 0, 0, 2.0 * x * factor - 1.0, false) *
                   std::pow(1.0 - y, i);

      const double y_contribution =
        Polynomials::jacobi_polynomial_value<double>(j, 2 * i + 1, 0, y, true);

      const double phi = x_contribution * y_contribution;

      if (std::fabs(phi) < 1e-14)
        return 0.0;

      return phi;
    }
  else if constexpr (dim == 3)
    {
      const double y = p[1];
      const double z = p[2];

      const double factor_yz =
        std::abs(1.0 - y - z) < 1e-14 ? 1.0 : 1.0 / (1.0 - y - z);
      const double factor_z = std::abs(1.0 - z) < 1e-14 ? 1.0 : 1.0 / (1.0 - z);

      const double x_contribution =
        i == 0 ? 1.0 :
                 Polynomials::jacobi_polynomial_value<double>(
                   i, 0, 0, 2.0 * x * factor_yz - 1.0, false) *
                   std::pow(1.0 - y - z, i);

      const double y_contribution =
        j == 0 ? 1.0 :
                 Polynomials::jacobi_polynomial_value<double>(
                   j, 2 * i + 1, 0, 2.0 * y * factor_z - 1.0, false) *
                   std::pow(1.0 - z, j);

      const double z_contribution =
        Polynomials::jacobi_polynomial_value<double>(
          k, 2 * (i + j) + 2, 0, z, true);

      const double phi = x_contribution * y_contribution * z_contribution;

      if (std::fabs(phi) < 1e-14)
        return 0.0;

      return phi;
    }

  DEAL_II_ASSERT_UNREACHABLE();
  return 0;
}



template <int dim>
double
ScalarLagrangePolynomialSimplex<dim>::evaluate_orthogonal_basis_function(
  const unsigned int i,
  const Point<dim>  &p) const
{
  AssertIndexRange(i, this->n());

  if constexpr (dim == 1)
    {
      // return the value of the orthogonal basis
      return evaluate_orthogonal_basis_function_by_degree(i, 0, 0, p);
    }
  else if constexpr (dim == 2)
    {
      // find corresponding entry to i
      // it holds 0 <= j + k <= degree
      for (unsigned int j = 0, counter = 0; j < this->degree() + 1; ++j)
        for (unsigned int k = 0; k < this->degree() + 1 - j; ++k, ++counter)
          if (counter == i)
            return evaluate_orthogonal_basis_function_by_degree(j, k, 0, p);
    }
  else if constexpr (dim == 3)
    {
      // find corresponding entry to i
      // it holds 0 <= j + k + l <= degree
      for (unsigned int j = 0, counter = 0; j < this->degree() + 1; ++j)
        for (unsigned int k = 0; k < this->degree() + 1 - j; ++k)
          for (unsigned int l = 0; l < this->degree() + 1 - j - k;
               ++l, ++counter)
            if (counter == i)
              return evaluate_orthogonal_basis_function_by_degree(j, k, l, p);
    }

  DEAL_II_ASSERT_UNREACHABLE();
  return 0;
}



template <int dim>
Tensor<1, dim>
ScalarLagrangePolynomialSimplex<dim>::
  evaluate_orthogonal_basis_derivative_by_degree(const unsigned int i,
                                                 const unsigned int j,
                                                 const unsigned int k,
                                                 const Point<dim>  &p) const
{
  AssertIndexRange(i + j + k, this->degree() + 1);

  Tensor<1, dim> grad;

  const double x = p[0];

  if constexpr (dim == 1)
    {
      grad[0] = Polynomials::jacobi_polynomial_derivative<double>(
                  i, 0, 0, 2.0 * x - 1.0, false) *
                2.0;
    }
  else if constexpr (dim == 2)
    {
      const double y = p[1];

      const double factor = std::abs(1.0 - y) < 1e-14 ? 0.0 : 1.0 / (1.0 - y);

      const double x_contribution =
        i == 0 ? 1.0 :
                 Polynomials::jacobi_polynomial_value<double>(
                   i, 0, 0, 2.0 * x * factor - 1.0, false) *
                   std::pow(1.0 - y, i);

      const double y_contribution =
        Polynomials::jacobi_polynomial_value<double>(j, 2 * i + 1, 0, y, true);


      const double x_derivative =
        i == 1 ? 2.0 :
                 Polynomials::jacobi_polynomial_derivative<double>(
                   i, 0, 0, 2.0 * x * factor - 1.0, false) *
                   std::pow(1.0 - y, std::max(1U, i)) * 2.0 * factor;

      grad[0] = x_derivative * y_contribution;

      const double y_derivative_x1 =
        i == 1 ?
          1.0 :
          Polynomials::jacobi_polynomial_derivative<double>(
            i, 0, 0, 2.0 * x * factor - 1.0, false) *
            2.0 * x * factor * factor * std::pow(1.0 - y, std::max(1U, i));

      const double y_derivative_x2 =
        i == 1 ? 0.0 :
                 Polynomials::jacobi_polynomial_value<double>(
                   i, 0, 0, 2.0 * x * factor - 1.0, false) *
                   i * (-1.0) * std::pow(1.0 - y, std::max(i - 1, 1U));

      const double y_derivative =
        Polynomials::jacobi_polynomial_derivative<double>(
          j, 2 * i + 1, 0, 2.0 * y - 1.0, false) *
        2.0;


      grad[1] = y_derivative_x1 * y_contribution +
                y_derivative_x2 * y_contribution +
                x_contribution * y_derivative;
    }
  else if constexpr (dim == 3)
    {
      const double y = p[1];
      const double z = p[2];

      const double factor_yz =
        std::abs(1.0 - y - z) < 1e-14 ? 0.0 : 1.0 / (1.0 - y - z);
      const double factor_z = std::abs(1.0 - z) < 1e-14 ? 0.0 : 1.0 / (1.0 - z);

      const double x_contribution =
        i == 0 ? 1.0 :
                 Polynomials::jacobi_polynomial_value<double>(
                   i, 0, 0, 2.0 * x * factor_yz - 1.0, false) *
                   std::pow(1.0 - y - z, i);

      const double y_contribution =
        j == 0 ? 1.0 :
                 Polynomials::jacobi_polynomial_value<double>(
                   j, 2 * i + 1, 0, 2.0 * y * factor_z - 1.0, false) *
                   std::pow(1.0 - z, j);

      const double z_contribution =
        Polynomials::jacobi_polynomial_value<double>(
          k, 2 * (i + j) + 2, 0, z, true);


      const double x_derivative =
        i == 1 ? 2.0 :
                 Polynomials::jacobi_polynomial_derivative<double>(
                   i, 0, 0, 2.0 * x * factor_yz - 1.0, false) *
                   std::pow(1.0 - y - z, std::max(1U, i)) * 2.0 * factor_yz;

      grad[0] = x_derivative * y_contribution * z_contribution;

      const double y_derivative_x1 =
        i == 1 ? 1.0 :
                 Polynomials::jacobi_polynomial_derivative<double>(
                   i, 0, 0, 2.0 * x * factor_yz - 1.0, false) *
                   2.0 * x * factor_yz * factor_yz *
                   std::pow(1.0 - y - z, std::max(1U, i));

      const double y_derivative_x2 =
        i == 1 ? 0.0 :
                 Polynomials::jacobi_polynomial_value<double>(
                   i, 0, 0, 2.0 * x * factor_yz - 1.0, false) *
                   i * (-1.0) * std::pow(1.0 - y - z, std::max(i - 1, 1U));

      const double y_derivative =
        j == 1 ? 2.0 * i + 3.0 :
                 Polynomials::jacobi_polynomial_derivative<double>(
                   j, 2 * i + 1, 0, 2.0 * y * factor_z - 1.0, false) *
                   std::pow(1.0 - z, std::max(1U, j)) * 2.0 * factor_z;

      grad[1] = y_derivative_x1 * y_contribution * z_contribution +
                y_derivative_x2 * y_contribution * z_contribution +
                x_contribution * y_derivative * z_contribution;


      const double z_derivative_y1 =
        j == 1 ?
          0.5 * (2.0 * i + 3.0) - 0.5 * (2.0 * i + 1.0) :
          Polynomials::jacobi_polynomial_derivative<double>(
            j, 2 * i + 1, 0, 2.0 * y * factor_z - 1.0, false) *
            2.0 * y * factor_z * factor_z * std::pow(1.0 - z, std::max(j, 1U));

      const double z_derivative_y2 =
        j == 1 ? 0.0 :
                 Polynomials::jacobi_polynomial_value<double>(
                   j, 2 * i + 1, 0, 2.0 * y * factor_z - 1.0, false) *
                   j * (-1.0) * std::pow(1.0 - z, std::max(1U, j - 1));

      const double z_derivative =
        Polynomials::jacobi_polynomial_derivative<double>(
          k, 2 * (i + j) + 2, 0, 2.0 * z - 1.0, false) *
        2.0;

      grad[2] = y_derivative_x1 * y_contribution * z_contribution +
                y_derivative_x2 * y_contribution * z_contribution +
                x_contribution * z_derivative_y1 * z_contribution +
                x_contribution * z_derivative_y2 * z_contribution +
                x_contribution * y_contribution * z_derivative;
    }
  else
    DEAL_II_NOT_IMPLEMENTED();

  for (unsigned int d = 0; d < dim; ++d)
    if (std::fabs(grad[d]) < 1e-14)
      grad[d] = 0.0;

  return grad;
}



template <int dim>
Tensor<1, dim>
ScalarLagrangePolynomialSimplex<dim>::evaluate_orthogonal_basis_derivative(
  const unsigned int i,
  const Point<dim>  &p) const
{
  AssertIndexRange(i, this->n());

  if constexpr (dim == 1)
    {
      // entrance to i corresponse to degree
      return evaluate_orthogonal_basis_derivative_by_degree(i, 0, 0, p);
    }
  else if constexpr (dim == 2)
    {
      // find corresponding entrance to i
      // it holds 0 <= j + k <= degree
      for (unsigned int j = 0, counter = 0; j < this->degree() + 1; ++j)
        for (unsigned int k = 0; k < this->degree() + 1 - j; ++k, ++counter)
          if (counter == i)
            return evaluate_orthogonal_basis_derivative_by_degree(j, k, 0, p);
    }
  else if constexpr (dim == 3)
    {
      // find corresponding entrance to i
      // it holds 0 <= j + k + l <= degree
      for (unsigned int j = 0, counter = 0; j < this->degree() + 1; ++j)
        for (unsigned int k = 0; k < this->degree() + 1 - j; ++k)
          for (unsigned int l = 0; l < this->degree() + 1 - j - k;
               ++l, ++counter)
            if (counter == i)
              return evaluate_orthogonal_basis_derivative_by_degree(j, k, l, p);
    }

  DEAL_II_ASSERT_UNREACHABLE();
  return Tensor<1, dim>();
}



template <int dim>
std::string
ScalarLagrangePolynomialSimplex<dim>::name() const
{
  return "ScalarLagrangePolynomialSimplex";
}



template <int dim>
std::unique_ptr<ScalarPolynomialsBase<dim>>
ScalarLagrangePolynomialSimplex<dim>::clone() const
{
  return std::make_unique<ScalarLagrangePolynomialSimplex<dim>>(*this);
}



template class ScalarLagrangePolynomialSimplex<1>;
template class ScalarLagrangePolynomialSimplex<2>;
template class ScalarLagrangePolynomialSimplex<3>;

DEAL_II_NAMESPACE_CLOSE
