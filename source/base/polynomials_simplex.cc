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
  // in 2D it holds that k = 0, again this multiplies by 1

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
Tensor<2, dim>
ScalarLagrangePolynomialSimplex<dim>::
  evaluate_orthogonal_basis_2nd_derivative_by_degree(const unsigned int i,
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

  Tensor<2, dim> deriv;

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
  // derivatives equal to zero

  // get the second derivatives over the product rule
  const double s     = 1 - z;
  const double ds_dz = -1.0;

  const double t     = 1 - y - z;
  const double dt_dy = -1.0;
  const double dt_dz = -1.0;

  // get the values of each polynomial
  const double Qi =
    Polynomials::jacobi_polynomial_homogenized_value<double>(i, 0, 0, x, t);
  const double Qj = Polynomials::jacobi_polynomial_homogenized_value<double>(
    j, 2 * i + 1, 0, y, s);
  const double Pk = Polynomials::jacobi_polynomial_value<double>(
    k, 2 * (i + j) + 2, 0, z, true);

  // get the first derivatives of Qi
  const double dQi_dx =
    Polynomials::jacobi_polynomial_homogenized_derivative<double>(
      1, 0, i, 0, 0, x, t);
  const double dQi_dt =
    Polynomials::jacobi_polynomial_homogenized_derivative<double>(
      0, 1, i, 0, 0, x, t);

  // get the first derivatives of Qj
  const double dQj_dy =
    Polynomials::jacobi_polynomial_homogenized_derivative<double>(
      1, 0, j, 2 * i + 1, 0, y, s);
  const double dQj_ds =
    Polynomials::jacobi_polynomial_homogenized_derivative<double>(
      0, 1, j, 2 * i + 1, 0, y, s);

  // get the first derivative of Pk
  const double dPk_dz = Polynomials::jacobi_polynomial_derivative<double>(
    k, 2 * (i + j) + 2, 0, z, true);

  // get the second and mixed derivatives of Qi
  const double dQi_dx_dx =
    Polynomials::jacobi_polynomial_homogenized_derivative<double>(
      2, 0, i, 0, 0, x, t);
  const double dQi_dx_dt =
    Polynomials::jacobi_polynomial_homogenized_derivative<double>(
      1, 1, i, 0, 0, x, t);
  const double dQi_dt_dt =
    Polynomials::jacobi_polynomial_homogenized_derivative<double>(
      0, 2, i, 0, 0, x, t);

  // get the second and mixed derivatives of Qj
  const double dQj_dy_dy =
    Polynomials::jacobi_polynomial_homogenized_derivative<double>(
      2, 0, j, 2 * i + 1, 0, y, s);
  const double dQj_dy_ds =
    Polynomials::jacobi_polynomial_homogenized_derivative<double>(
      1, 1, j, 2 * i + 1, 0, y, s);
  const double dQj_ds_ds =
    Polynomials::jacobi_polynomial_homogenized_derivative<double>(
      0, 2, j, 2 * i + 1, 0, y, s);

  // get the second derivative of Pk
  const double dPk_dz_dz =
    Polynomials::jacobi_polynomial_kth_derivative<double>(
      2, k, 2 * (i + j) + 2, 0, z, true);

  // now compute the entries
  deriv[0][0] = dQi_dx_dx * Qj * Pk;

  if constexpr (dim > 1)
    {
      deriv[0][1] = dQi_dx_dt * dt_dy * Qj * Pk + dQi_dx * dQj_dy * Pk;
      deriv[1][0] = deriv[0][1];
      deriv[1][1] = dQi_dt_dt * Qj * Pk + dQi_dt * dt_dy * dQj_dy * Pk * 2.0 +
                    Qi * dQj_dy_dy * Pk;

      if constexpr (dim > 2)
        {
          deriv[1][2] =
            dQi_dt_dt * Qj * Pk + dQi_dt * dt_dy * dQj_ds * ds_dz * Pk +
            dQi_dt * dt_dy * Qj * dPk_dz + dQi_dt * dt_dy * dQj_dy * Pk +
            Qi * dQj_dy_ds * ds_dz * Pk + Qi * dQj_dy * dPk_dz;
          deriv[2][1] = deriv[1][2];

          deriv[2][0] = dQi_dx_dt * dt_dz * Qj * Pk +
                        dQi_dx * dQj_ds * ds_dz * Pk + dQi_dx * Qj * dPk_dz;
          deriv[0][2] = deriv[2][0];

          deriv[2][2] =
            dQi_dt_dt * Qj * Pk + 2.0 * dQi_dt * dt_dz * dQj_ds * ds_dz * Pk +
            2.0 * dQi_dt * dt_dz * Qj * dPk_dz + Qi * dQj_ds_ds * Pk +
            2.0 * Qi * dQj_ds * ds_dz * dPk_dz + Qi * Qj * dPk_dz_dz;
        }
    }

  for (unsigned int d = 0; d < dim; ++d)
    for (unsigned int e = 0; e < dim; ++e)
      if (std::fabs(deriv[d][e]) < 1e-14)
        deriv[d][e] = 0.0;

  return deriv;
}



template <int dim>
Tensor<2, dim>
ScalarLagrangePolynomialSimplex<dim>::evaluate_orthogonal_basis_2nd_derivative(
  const unsigned int i,
  const Point<dim>  &p) const
{
  AssertIndexRange(i, this->n());

  if constexpr (dim == 1)
    {
      // entrance to i corresponse to degree
      return evaluate_orthogonal_basis_2nd_derivative_by_degree(i, 0, 0, p);
    }
  else if constexpr (dim == 2)
    {
      // find corresponding entrance to i
      // it holds 0 <= j + k <= degree
      for (unsigned int j = 0, counter = 0; j < this->degree() + 1; ++j)
        for (unsigned int k = 0; k < this->degree() + 1 - j; ++k, ++counter)
          if (counter == i)
            return evaluate_orthogonal_basis_2nd_derivative_by_degree(j,
                                                                      k,
                                                                      0,
                                                                      p);
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
              return evaluate_orthogonal_basis_2nd_derivative_by_degree(j,
                                                                        k,
                                                                        l,
                                                                        p);
    }

  DEAL_II_ASSERT_UNREACHABLE();
  return Tensor<2, dim>();
}



template <int dim>
Tensor<3, dim>
ScalarLagrangePolynomialSimplex<dim>::evaluate_orthogonal_basis_3rd_derivative(
  const unsigned int i,
  const Point<dim>  &p) const
{
  return evaluate_orthogonal_basis_nth_derivative<3>(i, p);
}



template <int dim>
Tensor<4, dim>
ScalarLagrangePolynomialSimplex<dim>::evaluate_orthogonal_basis_4th_derivative(
  const unsigned int i,
  const Point<dim>  &p) const
{
  return evaluate_orthogonal_basis_nth_derivative<4>(i, p);
}



template <int dim>
template <int derivative_order>
Tensor<derivative_order, dim>
ScalarLagrangePolynomialSimplex<dim>::evaluate_orthogonal_basis_nth_derivative(
  const unsigned int i,
  const Point<dim>  &p) const
{
  AssertIndexRange(i, this->n());

  if constexpr (dim == 1)
    {
      // entrance to i corresponse to degree
      return evaluate_orthogonal_basis_nth_derivative_by_degree<
        derivative_order>(i, 0, 0, p);
    }
  else if constexpr (dim == 2)
    {
      // find corresponding entrance to i
      // it holds 0 <= j + k <= degree
      for (unsigned int j = 0, counter = 0; j < this->degree() + 1; ++j)
        for (unsigned int k = 0; k < this->degree() + 1 - j; ++k, ++counter)
          if (counter == i)
            return evaluate_orthogonal_basis_nth_derivative_by_degree<
              derivative_order>(j, k, 0, p);
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
              return evaluate_orthogonal_basis_nth_derivative_by_degree<
                derivative_order>(j, k, l, p);
    }

  DEAL_II_ASSERT_UNREACHABLE();
  return Tensor<derivative_order, dim>();
}



template <int dim>
template <int derivative_order>
Tensor<derivative_order, dim>
ScalarLagrangePolynomialSimplex<dim>::
  evaluate_orthogonal_basis_nth_derivative_by_degree(const unsigned int i,
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

  Tensor<derivative_order, dim> deriv;

  // define t = 1 - y - z and s = 1 - z then
  // P_i^{0,0}(2x/t-1) * t^i
  // P_j^{2*i+1,0}(2*y/s-1)*s^j
  // P_k^{2 (i+j)+2,0}(2 z - 1)
  // =
  // Q_i^{0,0}(x,t)
  // Q_j^{2*i+1,0}(y,s)
  // P_k^{2*(i+j)+2,0}(2 z - 1)
  const double x = p[0];
  const double y = dim > 1 ? p[1] : 0.0;
  const double z = dim > 2 ? p[2] : 0.0;
  const double s = 1 - z;
  const double t = 1 - y - z;

  // helper function for small values of factorial
  // no need to get values higher than 4 as it is the highest derivative
  // calculated
  auto factorial = [](const unsigned int n) -> double {
    if (n <= 1)
      return 1.0;
    else if (n == 2)
      return 2.0;
    else if (n == 3)
      return 6.0;
    else if (n == 4)
      return 24.0;
    else
      DEAL_II_ASSERT_UNREACHABLE();
    return 0.0;
  };

  // helper function that computes the derivatives of the polynomials given the
  // indices of deriv
  auto compute_derivative_of_polynomials =
    [&](const std::array<unsigned int, derivative_order> indices) {
      // use Leibniz ruke to get the derivatives
      // https://en.wikipedia.org/wiki/General_Leibniz_rule
      // get the multi index of the derivative
      // e.g. deriv[0][0][0] is the third derivative with respect to x so
      // a_x = 3, a_y = 0, a_z = 0
      unsigned int a_x = 0;
      unsigned int a_y = 0;
      unsigned int a_z = 0;

      for (const auto idx : indices)
        {
          if (idx == 0)
            ++a_x;
          else if (idx == 1)
            ++a_y;
          else if (idx == 2)
            ++a_z;
          else
            DEAL_II_ASSERT_UNREACHABLE();
        }

      // needs to be derivative_order overall
      Assert(a_x + a_y + a_z == derivative_order, ExcInternalError());

      double derivative = 0.0;
      // now compute the derivative with multi index a
      // the multi index beta is for Qi
      // the multi index gamma is for Qj
      // the multi index delta is for Pk
      // sum over all multi indices, but some can be skipped
      // Qi is the only term dependent on x so beta_x = a_x
      // Qj is independent on x so gamma_x = 0
      // Pk is independent on x,y so delta_x, delta_y = 0
      for (unsigned int gamma_y = 0; gamma_y <= a_y; ++gamma_y)
        for (unsigned int gamma_z = 0; gamma_z <= a_z; ++gamma_z)
          for (unsigned int delta_z = 0; delta_z <= a_z - gamma_z; ++delta_z)
            {
              // from above
              const unsigned int beta_x = a_x;
              // derivative in y has to add up to a_y
              const unsigned int beta_y = a_y - gamma_y;
              // derivative in z has to add up to a_z
              const unsigned int beta_z = a_z - gamma_z - delta_z;

              const double coeff =
                factorial(a_y) * factorial(a_z) /
                (factorial(beta_y) * factorial(beta_z) * factorial(gamma_y) *
                 factorial(gamma_z) * factorial(delta_z));

              const double dQi_da =
                Polynomials::jacobi_polynomial_homogenized_derivative<double>(
                  beta_x, beta_y + beta_z, i, 0, 0, x, t) *
                std::pow(-1.0, beta_y + beta_z);

              const double dQj_da =
                Polynomials::jacobi_polynomial_homogenized_derivative<double>(
                  gamma_y, gamma_z, j, 2 * i + 1, 0, y, s) *
                std::pow(-1.0, gamma_z);

              const double dPk_da =
                delta_z == 0 ?
                  Polynomials::jacobi_polynomial_value<double>(
                    k, 2 * (i + j) + 2, 0, z, true) :
                  Polynomials::jacobi_polynomial_kth_derivative<double>(
                    delta_z, k, 2 * (i + j) + 2, 0, z, true);

              derivative += coeff * dQi_da * dQj_da * dPk_da;
            }

      if (std::abs(derivative) < 1e-12)
        return 0.0;

      return derivative;
    };

  if constexpr (derivative_order == 1)
    {
      // use the fast version here
      return evaluate_orthogonal_basis_derivative_by_degree(i, j, k, p);
    }
  else if constexpr (derivative_order == 2)
    {
      return evaluate_orthogonal_basis_2nd_derivative_by_degree(i, j, k, p);
    }
  else if constexpr (derivative_order == 3)
    {
      // loop over all entries in deriv
      for (unsigned int d1 = 0; d1 < dim; ++d1)
        for (unsigned int d2 = 0; d2 < dim; ++d2)
          for (unsigned int d3 = 0; d3 < dim; ++d3)
            {
              std::array<unsigned int, derivative_order> d = {{d1, d2, d3}};
              deriv[d1][d2][d3] = compute_derivative_of_polynomials(d);
            }
    }
  else if constexpr (derivative_order == 4)
    {
      // loop over all entries in deriv
      for (unsigned int d1 = 0; d1 < dim; ++d1)
        for (unsigned int d2 = 0; d2 < dim; ++d2)
          for (unsigned int d3 = 0; d3 < dim; ++d3)
            for (unsigned int d4 = 0; d4 < dim; ++d4)
              {
                std::array<unsigned int, derivative_order> d = {
                  {d1, d2, d3, d4}};
                deriv[d1][d2][d3][d4] = compute_derivative_of_polynomials(d);
              }
    }
  else
    DEAL_II_NOT_IMPLEMENTED();

  return deriv;
}



template <int dim>
void
ScalarLagrangePolynomialSimplex<dim>::evaluate(
  const Point<dim>            &unit_point,
  std::vector<double>         &values,
  std::vector<Tensor<1, dim>> &grads,
  std::vector<Tensor<2, dim>> &grad_grads,
  std::vector<Tensor<3, dim>> &third_derivatives,
  std::vector<Tensor<4, dim>> &fourth_derivatives) const
{
  if (values.size() == this->n())
    for (unsigned int i = 0; i < this->n(); ++i)
      values[i] = this->compute_value(i, unit_point);

  if (grads.size() == this->n())
    for (unsigned int i = 0; i < this->n(); ++i)
      grads[i] = this->compute_1st_derivative(i, unit_point);

  if (grad_grads.size() == this->n())
    for (unsigned int i = 0; i < this->n(); ++i)
      grad_grads[i] = this->compute_2nd_derivative(i, unit_point);

  if (third_derivatives.size() == this->n())
    for (unsigned int i = 0; i < this->n(); ++i)
      third_derivatives[i] = this->compute_3rd_derivative(i, unit_point);

  if (fourth_derivatives.size() == this->n())
    for (unsigned int i = 0; i < this->n(); ++i)
      fourth_derivatives[i] = this->compute_4th_derivative(i, unit_point);
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
