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
        else if (degree == 2)
          return 14;
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
  AssertIndexRange(this->degree(), 3);

  if (this->degree() == 1)
    {
      const double Q14 = 0.25;
      double       ration;

      const double r = p[0];
      const double s = p[1];
      const double t = p[2];

      if (fabs(t - 1.0) > 1.0e-14)
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
  else if (this->degree() == 2)
    {
      const double x = p[0];
      const double y = p[1];
      const double z = p[2];

      double result;

      double ration;
      double ration_square;

      if (fabs(z - 1.0) > 1.0e-14)
        {
          ration        = 1.0 / (z - 1.0);
          ration_square = std::pow(ration, 2);
        }
      else
        {
          ration        = 1.0;
          ration_square = 1.0;
        }

      switch (i)
        {
          case 0:
            result =
              (std::pow(z - 1, 2) *
                 ((1.0 / 12) * std::pow(x, 2) + (1.0 / 18) * x * (6 * z - 1) -
                  (1.0 / 36) * x + (1.0 / 12) * std::pow(y, 2) +
                  (1.0 / 18) * y * (6 * z - 1) - (1.0 / 36) * y +
                  (2.0 / 9) * std::pow(z, 2) - (7.0 / 36) * z - (1.0 / 36)) +
               (z - 1) * ((1.0 / 12) * x * y * (6 * z - 1) - (1.0 / 6) * x * y +
                          (1.0 / 12) * x *
                            (3 * std::pow(y, 2) - std::pow(z, 2) + 2 * z - 1) +
                          (1.0 / 12) * y *
                            (3 * std::pow(x, 2) - std::pow(z, 2) + 2 * z - 1)) +
               (1.0 / 36) * (3 * std::pow(x, 2) - std::pow(z, 2) + 2 * z - 1) *
                 (3 * std::pow(y, 2) - std::pow(z, 2) + 2 * z - 1)) *
              ration_square;
            break;
          case 1:
            result =
              (std::pow(z - 1, 2) *
                 ((1.0 / 12) * std::pow(x, 2) - (1.0 / 18) * x * (6 * z - 1) +
                  (1.0 / 36) * x + (1.0 / 12) * std::pow(y, 2) +
                  (1.0 / 18) * y * (6 * z - 1) - (1.0 / 36) * y +
                  (2.0 / 9) * std::pow(z, 2) - (7.0 / 36) * z - (1.0 / 36)) -
               (z - 1) * ((1.0 / 12) * x * y * (6 * z - 1) - (1.0 / 6) * x * y +
                          (1.0 / 12) * x *
                            (3 * std::pow(y, 2) - std::pow(z, 2) + 2 * z - 1) -
                          (1.0 / 12) * y *
                            (3 * std::pow(x, 2) - std::pow(z, 2) + 2 * z - 1)) +
               (1.0 / 36) * (3 * std::pow(x, 2) - std::pow(z, 2) + 2 * z - 1) *
                 (3 * std::pow(y, 2) - std::pow(z, 2) + 2 * z - 1)) *
              ration_square;
            break;
          case 2:
            result =
              (std::pow(z - 1, 2) *
                 ((1.0 / 12) * std::pow(x, 2) + (1.0 / 18) * x * (6 * z - 1) -
                  (1.0 / 36) * x + (1.0 / 12) * std::pow(y, 2) -
                  (1.0 / 18) * y * (6 * z - 1) + (1.0 / 36) * y +
                  (2.0 / 9) * std::pow(z, 2) - (7.0 / 36) * z - (1.0 / 36)) -
               (z - 1) * ((1.0 / 12) * x * y * (6 * z - 1) - (1.0 / 6) * x * y -
                          (1.0 / 12) * x *
                            (3 * std::pow(y, 2) - std::pow(z, 2) + 2 * z - 1) +
                          (1.0 / 12) * y *
                            (3 * std::pow(x, 2) - std::pow(z, 2) + 2 * z - 1)) +
               (1.0 / 36) * (3 * std::pow(x, 2) - std::pow(z, 2) + 2 * z - 1) *
                 (3 * std::pow(y, 2) - std::pow(z, 2) + 2 * z - 1)) *
              ration_square;
            break;
          case 3:
            result =
              (std::pow(z - 1, 2) *
                 ((1.0 / 12) * std::pow(x, 2) - (1.0 / 18) * x * (6 * z - 1) +
                  (1.0 / 36) * x + (1.0 / 12) * std::pow(y, 2) -
                  (1.0 / 18) * y * (6 * z - 1) + (1.0 / 36) * y +
                  (2.0 / 9) * std::pow(z, 2) - (7.0 / 36) * z - (1.0 / 36)) -
               (z - 1) *
                 (-(1.0 / 12) * x * y * (6 * z - 1) + (1.0 / 6) * x * y +
                  (1.0 / 12) * x *
                    (3 * std::pow(y, 2) - std::pow(z, 2) + 2 * z - 1) +
                  (1.0 / 12) * y *
                    (3 * std::pow(x, 2) - std::pow(z, 2) + 2 * z - 1)) +
               (1.0 / 36) * (3 * std::pow(x, 2) - std::pow(z, 2) + 2 * z - 1) *
                 (3 * std::pow(y, 2) - std::pow(z, 2) + 2 * z - 1)) *
              ration_square;
            break;
          case 4:
            result = z * (2 * z - 1);
            break;
          case 5:
            result =
              (-x * (z - 1) *
                 (0.5 * std::pow(y, 2) - 1.0 / 6 * std::pow(z, 2) +
                  1.0 / 3 * z - 1.0 / 6) +
               std::pow(z - 1, 2) *
                 (1.0 / 3 * std::pow(x, 2) + 1.0 / 18 * x * (6 * z - 1) -
                  5.0 / 18 * x - 1.0 / 6 * std::pow(y, 2) +
                  1.0 / 18 * std::pow(z, 2) - 1.0 / 9 * z + 1.0 / 18) -
               1.0 / 18 * (3 * std::pow(x, 2) - std::pow(z, 2) + 2 * z - 1) *
                 (3 * std::pow(y, 2) - std::pow(z, 2) + 2 * z - 1)) *
              ration_square;
            break;
          case 6:
            result =
              (x * (z - 1) *
                 (0.5 * std::pow(y, 2) - 1.0 / 6 * std::pow(z, 2) +
                  1.0 / 3 * z - 1.0 / 6) +
               std::pow(z - 1, 2) *
                 (1.0 / 3 * std::pow(x, 2) - 1.0 / 18 * x * (6 * z - 1) +
                  5.0 / 18 * x - 1.0 / 6 * std::pow(y, 2) +
                  1.0 / 18 * std::pow(z, 2) - 1.0 / 9 * z + 1.0 / 18) -
               1.0 / 18 * (3 * std::pow(x, 2) - std::pow(z, 2) + 2 * z - 1) *
                 (3 * std::pow(y, 2) - std::pow(z, 2) + 2 * z - 1)) *
              ration_square;
            break;
          case 7:
            result =
              (-1.0 / 6 * y * (z - 1) *
                 (3 * std::pow(x, 2) - std::pow(z, 2) + 2 * z - 1) +
               std::pow(z - 1, 2) *
                 (-1.0 / 6 * std::pow(x, 2) + 1.0 / 3 * std::pow(y, 2) +
                  1.0 / 18 * y * (6 * z - 1) - 5.0 / 18 * y +
                  1.0 / 18 * std::pow(z, 2) - 1.0 / 9 * z + 1.0 / 18) -
               1.0 / 18 * (3 * std::pow(x, 2) - std::pow(z, 2) + 2 * z - 1) *
                 (3 * std::pow(y, 2) - std::pow(z, 2) + 2 * z - 1)) *
              ration_square;
            break;
          case 8:
            result =
              (1.0 / 6 * y * (z - 1) *
                 (3 * std::pow(x, 2) - std::pow(z, 2) + 2 * z - 1) +
               std::pow(z - 1, 2) *
                 (-1.0 / 6 * std::pow(x, 2) + 1.0 / 3 * std::pow(y, 2) -
                  1.0 / 18 * y * (6 * z - 1) + 5.0 / 18 * y +
                  1.0 / 18 * std::pow(z, 2) - 1.0 / 9 * z + 1.0 / 18) -
               1.0 / 18 * (3 * std::pow(x, 2) - std::pow(z, 2) + 2 * z - 1) *
                 (3 * std::pow(y, 2) - std::pow(z, 2) + 2 * z - 1)) *
              ration_square;
            break;
          case 9:
            result =
              (-x * y * z - x * std::pow(z, 2) + x * z - y * std::pow(z, 2) +
               y * z - std::pow(z, 3) + 2.0 * std::pow(z, 2) - z - 0) *
              ration;
            break;
          case 10:
            result =
              (x * y * z + x * std::pow(z, 2) - x * z - y * std::pow(z, 2) +
               y * z - std::pow(z, 3) + 2.0 * std::pow(z, 2) - z - 0) *
              ration;
            break;
          case 11:
            result =
              (x * y * z - x * std::pow(z, 2) + x * z + y * std::pow(z, 2) -
               y * z - std::pow(z, 3) + 2.0 * std::pow(z, 2) - z + 0) *
              ration;
            break;
          case 12:
            result = z *
                     (-x * y + x * z - x + y * z - y - std::pow(z, 2) +
                      2.0 * z - 1.0) *
                     ration;
            break;
          case 13:
            result =
              (std::pow(z - 1, 2) *
                 (-2.0 / 3 * std::pow(x, 2) - 2.0 / 3 * std::pow(y, 2) +
                  8.0 / 9 * std::pow(z, 2) - 16.0 / 9 * z + 8.0 / 9) +
               1.0 / 9 * (3 * std::pow(x, 2) - std::pow(z, 2) + 2 * z - 1) *
                 (3 * std::pow(y, 2) - std::pow(z, 2) + 2 * z - 1)) *
              ration_square;
            break;
          default:
            DEAL_II_ASSERT_UNREACHABLE();
        }
      return result;
    }
  DEAL_II_ASSERT_UNREACHABLE();
  return 0.0;
}



template <int dim>
Tensor<1, dim>
ScalarLagrangePolynomialPyramid<dim>::compute_grad(const unsigned int i,
                                                   const Point<dim>  &p) const
{
  AssertDimension(dim, 3);
  AssertIndexRange(this->degree(), 3);

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

      if (fabs(t - 1.0) > 1.0e-14)
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
  else if (this->degree() == 2)
    {
      const double x = p[0];
      const double y = p[1];
      const double z = p[2];

      switch (i)
        {
          case 0:
            grad[0] = (1.0 / 6.0 * x *
                         (3 * std::pow(y, 2) - std::pow(z, 2) + 2 * z - 1) +
                       std::pow((z - 1), 2) *
                         (1.0 / 6.0 * x + 1.0 / 3.0 * z - 1.0 / 12.0) +
                       (z - 1) * (0.5 * x * y + 0.25 * std::pow(y, 2) +
                                  1.0 / 12.0 * y * (6 * z - 1) - 1.0 / 6.0 * y -
                                  1.0 / 12.0 * std::pow(z, 2) + 1.0 / 6.0 * z -
                                  1.0 / 12.0)) /
                      std::pow((z - 1), 2);
            grad[1] = (1.0 / 6.0 * y *
                         (3 * std::pow(x, 2) - std::pow(z, 2) + 2 * z - 1) +
                       std::pow((z - 1), 2) *
                         (1.0 / 6.0 * y + 1.0 / 3.0 * z - 1.0 / 12.0) +
                       (z - 1) * (0.25 * std::pow(x, 2) + 0.5 * x * y +
                                  1.0 / 12.0 * x * (6 * z - 1) - 1.0 / 6.0 * x -
                                  1.0 / 12.0 * std::pow(z, 2) + 1.0 / 6.0 * z -
                                  1.0 / 12.0)) /
                      std::pow((z - 1), 2);
            grad[2] =
              (-0.5 * std::pow(x, 2) * std::pow(y, 2) -
               0.25 * std::pow(x, 2) * y * z + 0.25 * std::pow(x, 2) * y -
               0.25 * x * std::pow(y, 2) * z + 0.25 * x * std::pow(y, 2) -
               0.25 * x * y * z + 0.25 * x * y + 0.25 * x * std::pow(z, 3) -
               0.75 * x * std::pow(z, 2) + 0.75 * x * z - 0.25 * x +
               0.25 * y * std::pow(z, 3) - 0.75 * y * std::pow(z, 2) +
               0.75 * y * z - 0.25 * y + 0.5 * std::pow(z, 4) -
               1.75 * std::pow(z, 3) + 2.25 * std::pow(z, 2) - 1.25 * z +
               0.25) /
              (std::pow(z, 3) - 3.0 * std::pow(z, 2) + 3.0 * z - 1.0);
            break;

          case 1:
            grad[0] = (1.0 / 6.0 * x *
                         (3 * std::pow(y, 2) - std::pow(z, 2) + 2 * z - 1) +
                       std::pow((z - 1), 2) *
                         (1.0 / 6.0 * x - 1.0 / 3.0 * z + 1.0 / 12.0) +
                       (z - 1) * (0.5 * x * y - 0.25 * std::pow(y, 2) -
                                  1.0 / 12.0 * y * (6 * z - 1) + 1.0 / 6.0 * y +
                                  1.0 / 12.0 * std::pow(z, 2) - 1.0 / 6.0 * z +
                                  1.0 / 12.0)) /
                      std::pow((z - 1), 2);
            grad[1] = (1.0 / 6.0 * y *
                         (3 * std::pow(x, 2) - std::pow(z, 2) + 2 * z - 1) +
                       std::pow((z - 1), 2) *
                         (1.0 / 6.0 * y + 1.0 / 3.0 * z - 1.0 / 12.0) -
                       (z - 1) * (-0.25 * std::pow(x, 2) + 0.5 * x * y -
                                  1.0 / 12.0 * x * (6 * z - 1) + 1.0 / 6.0 * x +
                                  1.0 / 12.0 * std::pow(z, 2) - 1.0 / 6.0 * z +
                                  1.0 / 12.0)) /
                      std::pow((z - 1), 2);
            grad[2] =
              (-0.5 * std::pow(x, 2) * std::pow(y, 2) -
               0.25 * std::pow(x, 2) * y * z + 0.25 * std::pow(x, 2) * y +
               0.25 * x * std::pow(y, 2) * z - 0.25 * x * std::pow(y, 2) +
               0.25 * x * y * z - 0.25 * x * y - 0.25 * x * std::pow(z, 3) +
               0.75 * x * std::pow(z, 2) - 0.75 * x * z + 0.25 * x +
               0.25 * y * std::pow(z, 3) - 0.75 * y * std::pow(z, 2) +
               0.75 * y * z - 0.25 * y + 0.5 * std::pow(z, 4) -
               1.75 * std::pow(z, 3) + 2.25 * std::pow(z, 2) - 1.25 * z +
               0.25) /
              (std::pow(z, 3) - 3.0 * std::pow(z, 2) + 3.0 * z - 1.0);
            break;

          case 2:
            grad[0] = (1.0 / 6.0 * x *
                         (3 * std::pow(y, 2) - std::pow(z, 2) + 2 * z - 1) +
                       std::pow((z - 1), 2) *
                         (1.0 / 6.0 * x + 1.0 / 3.0 * z - 1.0 / 12.0) -
                       (z - 1) * (0.5 * x * y - 0.25 * std::pow(y, 2) +
                                  1.0 / 12.0 * y * (6 * z - 1) - 1.0 / 6.0 * y +
                                  1.0 / 12.0 * std::pow(z, 2) - 1.0 / 6.0 * z +
                                  1.0 / 12.0)) /
                      std::pow((z - 1), 2);
            grad[1] = (1.0 / 6.0 * y *
                         (3 * std::pow(x, 2) - std::pow(z, 2) + 2 * z - 1) +
                       std::pow((z - 1), 2) *
                         (1.0 / 6.0 * y - 1.0 / 3.0 * z + 1.0 / 12.0) +
                       (z - 1) * (-0.25 * std::pow(x, 2) + 0.5 * x * y -
                                  1.0 / 12.0 * x * (6 * z - 1) + 1.0 / 6.0 * x +
                                  1.0 / 12.0 * std::pow(z, 2) - 1.0 / 6.0 * z +
                                  1.0 / 12.0)) /
                      std::pow((z - 1), 2);
            grad[2] =
              (-0.5 * std::pow(x, 2) * std::pow(y, 2) +
               0.25 * std::pow(x, 2) * y * z - 0.25 * std::pow(x, 2) * y -
               0.25 * x * std::pow(y, 2) * z + 0.25 * x * std::pow(y, 2) +
               0.25 * x * y * z - 0.25 * x * y + 0.25 * x * std::pow(z, 3) -
               0.75 * x * std::pow(z, 2) + 0.75 * x * z - 0.25 * x -
               0.25 * y * std::pow(z, 3) + 0.75 * y * std::pow(z, 2) -
               0.75 * y * z + 0.25 * y + 0.5 * std::pow(z, 4) -
               1.75 * std::pow(z, 3) + 2.25 * std::pow(z, 2) - 1.25 * z +
               0.25) /
              (std::pow(z, 3) - 3 * std::pow(z, 2) + 3 * z - 1.0);
            break;
          case 3:
            grad[0] = (1.0 / 6.0 * x *
                         (3 * std::pow(y, 2) - std::pow(z, 2) + 2 * z - 1) +
                       std::pow((z - 1), 2) *
                         (1.0 / 6.0 * x - 1.0 / 3.0 * z + 1.0 / 12.0) -
                       (z - 1) * (0.5 * x * y + 0.25 * std::pow(y, 2) -
                                  1.0 / 12.0 * y * (6 * z - 1) + 1.0 / 6.0 * y -
                                  1.0 / 12.0 * std::pow(z, 2) + 1.0 / 6.0 * z -
                                  1.0 / 12.0)) /
                      std::pow((z - 1), 2);
            grad[1] = (1.0 / 6.0 * y *
                         (3 * std::pow(x, 2) - std::pow(z, 2) + 2 * z - 1) +
                       std::pow((z - 1), 2) *
                         (1.0 / 6.0 * y - 1.0 / 3.0 * z + 1.0 / 12.0) -
                       (z - 1) * (0.25 * std::pow(x, 2) + 0.5 * x * y -
                                  1.0 / 12.0 * x * (6 * z - 1) + 1.0 / 6.0 * x -
                                  1.0 / 12.0 * std::pow(z, 2) + 1.0 / 6.0 * z -
                                  1.0 / 12.0)) /
                      std::pow((z - 1), 2);
            grad[2] =
              (-0.5 * std::pow(x, 2) * std::pow(y, 2) +
               0.25 * std::pow(x, 2) * y * z - 0.25 * std::pow(x, 2) * y +
               0.25 * x * std::pow(y, 2) * z - 0.25 * x * std::pow(y, 2) -
               0.25 * x * y * z + 0.25 * x * y - 0.25 * x * std::pow(z, 3) +
               0.75 * x * std::pow(z, 2) - 0.75 * x * z + 0.25 * x -
               0.25 * y * std::pow(z, 3) + 0.75 * y * std::pow(z, 2) -
               0.75 * y * z + 0.25 * y + 0.5 * std::pow(z, 4) -
               1.75 * std::pow(z, 3) + 2.25 * std::pow(z, 2) - 1.25 * z +
               0.25) /
              (std::pow(z, 3) - 3 * std::pow(z, 2) + 3 * z - 1.0);
            break;
          case 4:
            grad[0] = 0;
            grad[1] = 0;
            grad[2] =
              (4.0 * std::pow(z, 3) - 9.0 * std::pow(z, 2) + 6.0 * z - 1.0) /
              (std::pow(z, 2) - 2.0 * z + 1.0);
            break;
          case 5:
            grad[0] =
              (-1.0 / 3.0 * x *
                 (3 * std::pow(y, 2) - std::pow(z, 2) + 2 * z - 1) +
               std::pow((z - 1), 2) *
                 (2.0 / 3.0 * x + 1.0 / 3.0 * z - 1.0 / 3.0) -
               (z - 1) * (0.5 * std::pow(y, 2) - 1.0 / 6.0 * std::pow(z, 2) +
                          1.0 / 3.0 * z - 1.0 / 6.0)) /
              std::pow((z - 1), 2);
            grad[1] = (-x * (y * (z - 1)) -
                       1.0 / 3.0 * y *
                         (3 * std::pow(x, 2) - std::pow(z, 2) + 2 * z - 1) +
                       (-1.0 / 3.0 * y) * (std::pow((z - 1), 2))) /
                      std::pow((z - 1), 2);
            grad[2] =
              (std::pow(x, 2) * std::pow(y, 2) + 0.5 * x * std::pow(y, 2) * z -
               0.5 * x * std::pow(y, 2) + 0.5 * x * std::pow(z, 3) -
               1.5 * x * std::pow(z, 2) + 1.5 * x * z - 0.5 * x) /
              (std::pow(z, 3) - 3 * std::pow(z, 2) + 3 * z - 1.0);
            break;
          case 6:
            grad[0] = (-1.0 / 3.0 * x *
                         (3 * std::pow(y, 2) - std::pow(z, 2) + 2 * z - 1) +
                       std::pow((z - 1), 2) *
                         (2.0 / 3.0 * x - 1.0 / 3.0 * z + 1.0 / 3.0) +
                       (z - 1) * (0.5 * std::pow(y, 2) + 0 -
                                  1.0 / 6.0 * std::pow(z, 2) + 1.0 / 3.0 * z -
                                  1.0 / 6.0)) /
                      std::pow((z - 1), 2);
            grad[1] = (x * (1.0 * y + 0) * (z - 1) -
                       1.0 / 3.0 * y *
                         (3 * std::pow(x, 2) - std::pow(z, 2) + 2 * z - 1) -
                       (1.0 / 3.0 * y + 0) * std::pow((z - 1), 2)) /
                      std::pow((z - 1), 2);
            grad[2] =
              (1.0 * std::pow(x, 2) * std::pow(y, 2) + 0 - 0 + 0 -
               0.5 * x * std::pow(y, 2) * z + 0.5 * x * std::pow(y, 2) - 0 + 0 -
               0.5 * x * std::pow(z, 3) + 1.5 * x * std::pow(z, 2) -
               1.5 * x * z + 0.5 * x - 0 + 0 - 0 + 0 + 0 - 0 + 0) /
              (std::pow(z, 3) - 3.0 * std::pow(z, 2) + 3.0 * z - 1.0);
            break;
          case 7:
            grad[0] = x *
                      (-1.0 * std::pow(y, 2) - 1.0 * y * (z - 1) +
                       1.0 / 3.0 * std::pow(z, 2) - 2.0 / 3.0 * z -
                       1.0 / 3.0 * std::pow((z - 1), 2) + 1.0 / 3.0) /
                      std::pow((z - 1), 2);
            grad[1] = (-1.0 / 3.0 * y *
                         (3 * std::pow(x, 2) - std::pow(z, 2) + 2 * z - 1) +
                       std::pow((z - 1), 2) *
                         (2.0 / 3.0 * y + 1.0 / 3.0 * z - 1.0 / 3.0) -
                       1.0 / 6.0 * (z - 1) *
                         (3 * std::pow(x, 2) - std::pow(z, 2) + 2 * z - 1)) /
                      std::pow((z - 1), 2);
            grad[2] =
              (1.0 * std::pow(x, 2) * std::pow(y, 2) +
               0.5 * std::pow(x, 2) * y * z - 0.5 * std::pow(x, 2) * y + 0 - 0 +
               0 + 0.5 * y * std::pow(z, 3) - 1.5 * y * std::pow(z, 2) +
               1.5 * y * z - 0.5 * y - 0 - 0 + 0) /
              (std::pow(z, 3) - 3.0 * std::pow(z, 2) + 3.0 * z - 1.0);
            break;
          case 8:
            grad[0] = x *
                      (-1.0 * std::pow(y, 2) + 1.0 * y * (z - 1) +
                       1.0 / 3.0 * std::pow(z, 2) - 2.0 / 3.0 * z -
                       1.0 / 3.0 * std::pow((z - 1), 2) + 1.0 / 3.0) /
                      std::pow((z - 1), 2);
            grad[1] = (-1.0 / 3.0 * y *
                         (3 * std::pow(x, 2) - std::pow(z, 2) + 2 * z - 1) +
                       std::pow((z - 1), 2) *
                         (2.0 / 3.0 * y - 1.0 / 3.0 * z + 1.0 / 3.0) +
                       1.0 / 6.0 * (z - 1) *
                         (3 * std::pow(x, 2) - std::pow(z, 2) + 2 * z - 1)) /
                      std::pow((z - 1), 2);
            grad[2] = (1.0 * std::pow(x, 2) * std::pow(y, 2) -
                       0.5 * std::pow(x, 2) * y * z + 0.5 * std::pow(x, 2) * y +
                       0 - 0.5 * y * std::pow(z, 3) + 1.5 * y * std::pow(z, 2) -
                       1.5 * y * z + 0.5 * y - 0 - 0 + 0) /
                      (std::pow(z, 3) - 3.0 * std::pow(z, 2) + 3.0 * z - 1.0);
            break;
          case 9:
            grad[0] = z * (-y - z + 1) / (z - 1);
            grad[1] = z * (-x - z + 1) / (z - 1);
            grad[2] = (x * y - x * std::pow(z, 2) + 2 * x * z - x -
                       y * std::pow(z, 2) + 2 * y * z - y - 2 * std::pow(z, 3) +
                       5 * std::pow(z, 2) - 4 * z + 1) /
                      (std::pow(z, 2) - 2 * z + 1);
            break;
          case 10:
            grad[0] = z * (y + z - 1) / (z - 1);
            grad[1] = z * (x - z + 1) / (z - 1);
            grad[2] = (-x * y + x * std::pow(z, 2) - 2 * x * z + x -
                       y * std::pow(z, 2) + 2 * y * z - y - 2 * std::pow(z, 3) +
                       5 * std::pow(z, 2) - 4 * z + 1) /
                      (std::pow(z, 2) - 2 * z + 1);
            break;
          case 11:
            grad[0] = z * (y - z + 1) / (z - 1);
            grad[1] = z * (x + z - 1) / (z - 1);
            grad[2] = (-x * y - x * std::pow(z, 2) + 2 * x * z - x +
                       y * std::pow(z, 2) - 2 * y * z + y - 2 * std::pow(z, 3) +
                       5 * std::pow(z, 2) - 4 * z + 1) /
                      (std::pow(z, 2) - 2 * z + 1);
            break;
          case 12:
            grad[0] = z * (-y + z - 1) / (z - 1);
            grad[1] = z * (-x + z - 1) / (z - 1);
            grad[2] = (x * y + x * std::pow(z, 2) - 2 * x * z + x +
                       y * std::pow(z, 2) - 2 * y * z + y - 2 * std::pow(z, 3) +
                       5 * std::pow(z, 2) - 4 * z + 1) /
                      (std::pow(z, 2) - 2 * z + 1);
            break;
          case 13:
            grad[0] = x *
                      (2 * std::pow(y, 2) - 2 * std::pow(z, 2) + 4 * z - 2) /
                      (std::pow(z, 2) - 2 * z + 1);
            grad[1] =
              y *
              (2.0 * std::pow(x, 2) - 2.0 * std::pow(z, 2) + 4.0 * z - 2.0) /
              (std::pow(z, 2) - 2.0 * z + 1.0);
            grad[2] = (-2.0 * std::pow(x, 2) * std::pow(y, 2) - 0 + 0 - 0 - 0 +
                       0 + 2.0 * std::pow(z, 4) - 8.0 * std::pow(z, 3) +
                       12.0 * std::pow(z, 2) - 8.0 * z + 2.0) /
                      (std::pow(z, 3) - 3.0 * std::pow(z, 2) + 3.0 * z - 1.0);
            break;

          default:
            DEAL_II_ASSERT_UNREACHABLE();
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
