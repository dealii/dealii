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

  if constexpr (dim == 1)
    Assert(j == 0 && k == 0, ExcInternalError());
  else if constexpr (dim == 2)
    Assert(k == 0, ExcInternalError());
  else if constexpr (dim == 3)
    {
      // nothing to assert
    }
  else
    DEAL_II_ASSERT_UNREACHABLE();

  const double x = p[0];
  const double y = dim > 1 ? p[1] : 0.0;
  const double z = dim > 2 ? p[2] : 0.0;

  // the basis function looks like
  // P_i^{0,0}(2x/(1-y-z)-1) * (1-y-z)^i
  // P_j^{2*i+1,0}(2*y/(1-z)-1)*(1-z)^j
  // P_k^{2 (i+j)+2,0}(2 z - 1)
  // like in 2d use the homogenized Jacobi polynomials
  // define t = 1 - y - z and s = 1 - z
  // the first term becomes
  // Q_i^{0,0}(x,t) = P_i^{0,0}(2x/(1-y-z)-1) * (1-y-z)^i
  // and the second
  // Q_j^{2i+1,0}(y,s) = P_j^{2i+1,0}(2y/(1-z)-1) * (1-z)^j

  // in 1D it holds that j,k = 0, such the (homogenized) Jacobi polynomials
  // describing the y and z contributions just equal to 1
  // in 2D it hold that  k = 0, again this multiplies by 1

  const double s = 1 - z;
  const double t = 1 - y - z;

  const double Qi =
    Polynomials::jacobi_polynomial_homogenized_value<double>(i, 0, 0, x, t);

  const double Qj = Polynomials::jacobi_polynomial_homogenized_value<double>(
    j, 2 * i + 1, 0, y, s);

  const double Pk = Polynomials::jacobi_polynomial_value<double>(
    k, 2 * (i + j) + 2, 0, z, true);

  const double phi = Qi * Qj * Pk;

  if (std::fabs(phi) < 1e-14)
    return 0.0;

  return phi;
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

  if constexpr (dim == 1)
    Assert(j == 0 && k == 0, ExcInternalError());
  else if constexpr (dim == 2)
    Assert(k == 0, ExcInternalError());
  else if constexpr (dim == 3)
    {
      // nothing to assert
    }
  else
    DEAL_II_ASSERT_UNREACHABLE();

  Tensor<1, dim> grad;

  const double x = p[0];
  const double y = dim > 1 ? p[1] : 0.0;
  const double z = dim > 2 ? p[2] : 0.0;

  // define t = 1 - y - z and s = 1 - z then
  // P_i^{0,0}(2x/t-1) * t^i
  // P_j^{2*i+1,0}(2*y/s-1)*s^j
  // P_k^{2 (i+j)+2,0}(2 z - 1)
  // =
  // Q_i^{0,0}(x,t)
  // Q_j^{2*i+1,0}(y,s)
  // P_k^{2*(i+j)+2,0}(2 z - 1)

  // The 1D (2D) cases are again covered by having j,k = 0 (k = 0) such that
  // the contributions to the value equal to one and the contributions to the
  // derivative equal to zero

  // get the derivatives over the product rule
  const double s     = 1 - z;
  const double ds_dz = -1.0;

  const double t     = 1 - y - z;
  const double dt_dy = -1.0;
  const double dt_dz = -1.0;

  const double Qi =
    Polynomials::jacobi_polynomial_homogenized_value<double>(i, 0, 0, x, t);
  const double Qj = Polynomials::jacobi_polynomial_homogenized_value<double>(
    j, 2 * i + 1, 0, y, s);
  const double Pk = Polynomials::jacobi_polynomial_value<double>(
    k, 2 * (i + j) + 2, 0, z, true);

  const double dQi_dx =
    Polynomials::jacobi_polynomial_homogenized_derivative<double>(
      1, 0, i, 0, 0, x, t);
  const double dQi_dt =
    Polynomials::jacobi_polynomial_homogenized_derivative<double>(
      0, 1, i, 0, 0, x, t);

  const double dQj_dy =
    Polynomials::jacobi_polynomial_homogenized_derivative<double>(
      1, 0, j, 2 * i + 1, 0, y, s);
  const double dQj_ds =
    Polynomials::jacobi_polynomial_homogenized_derivative<double>(
      0, 1, j, 2 * i + 1, 0, y, s);

  const auto dPk_dz = Polynomials::jacobi_polynomial_derivative<double>(
    k, 2 * (i + j) + 2, 0, z, true);


  grad[0] = dQi_dx * Qj * Pk;
  if constexpr (dim > 1)
    grad[1] = dQi_dt * dt_dy * Qj * Pk + Qi * dQj_dy * Pk;
  if constexpr (dim > 2)
    grad[2] =
      dQi_dt * dt_dz * Qj * Pk + Qi * dQj_ds * ds_dz * Pk + Qi * Qj * dPk_dz;

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
