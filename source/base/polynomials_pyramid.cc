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
double
ScalarLagrangePolynomialPyramid<dim>::compute_jacobi_basis_functions(
  const unsigned int i,
  const Point<dim>  &p) const
{
  if (this->degree() == 2)
    {
      double       basis;
      const double x = p[0];
      const double y = p[1];
      const double z = p[2];

      double ratio;
      if (std::fabs(z - 1.0) < 1e-12)
        {
          ratio = 0.0;
        }
      else
        {
          ratio = 1.0 / (z - 1.0);
        }

      switch (i)
        {
          case 0:
            basis = 1.0;
            break;
          case 1:
            basis = 4.0 * z - 1.0;
            break;
          case 2:
            basis = 15.0 * std::pow(z, 2) - 10.0 * z + 1.0;
            break;
          case 3:
            basis = y;
            break;
          case 4:
            basis = y * (6.0 * z - 1.0);
            break;
          case 5:
            basis = 3 * std::pow(y, 2) * 0.5 - std::pow(z, 2) * 0.5 + z - 0.5;
            break;
          case 6:
            basis = x;
            break;
          case 7:
            basis = x * (6.0 * z - 1.0);
            break;
          case 8:
            basis = -x * y * ratio;
            break;
          case 9:
            basis = x * y * (1.0 - 6.0 * z) * ratio;
            break;
          case 10:
            basis = x * (-3.0 * std::pow(y, 2) + std::pow(z, 2) - 2 * z + 1.0) *
                    0.5 * ratio;
            break;
          case 11:
            basis = 3.0 / 2.0 * std::pow(x, 2) - std::pow(z, 2) * 0.5 + z - 0.5;
            break;
          case 12:
            basis = y *
                    (-3.0 * std::pow(x, 2) + std::pow(z, 2) - 2.0 * z + 1.0) *
                    0.5 * ratio;
            break;
          case 13:
            basis = (3.0 * std::pow(x, 2) - std::pow(z, 2) + 2.0 * z - 1.0) *
                    (3.0 * std::pow(y, 2) - std::pow(z, 2) + 2.0 * z - 1.0) *
                    0.25 * std::pow(ratio, 2);
            break;
          default:

            DEAL_II_NOT_IMPLEMENTED();
            break;
        }
      if (std::fabs(basis) < 1e-14)
        basis = 0.0;

      return basis;
    }
  DEAL_II_NOT_IMPLEMENTED();
  return 0;
}


template <int dim>
Tensor<1, dim>
ScalarLagrangePolynomialPyramid<dim>::compute_jacobi_deriv_basis_functions(
  const unsigned int i,
  const Point<dim>  &p) const
{
  if (this->degree() == 2)
    {
      Tensor<1, dim> grad;
      const double   x = p[0];
      const double   y = p[1];
      const double   z = p[2];

      double ratio;
      if (std::fabs(z - 1.0) < 1e-12)
        {
          ratio = 0.0;
        }
      else
        {
          ratio = 1.0 / (z - 1.0);
        }

      switch (i)
        {
          case 0:
            grad[0] = 0.0;
            grad[1] = 0.0;
            grad[2] = 0.0;
            break;
          case 1:
            grad[0] = 0.0;
            grad[1] = 0.0;
            grad[2] = 4.0;
            break;
          case 2:
            grad[0] = 0.0;
            grad[1] = 0.0;
            grad[2] = 30.0 * z - 10.0;
            break;
          case 3:
            grad[0] = 0.0;
            grad[1] = 1.0;
            grad[2] = 0.0;
            break;
          case 4:
            grad[0] = 0.0;
            grad[1] = 6.0 * z - 1.0;
            grad[2] = 6.0 * y;
            break;
          case 5:
            grad[0] = 0.0;
            grad[1] = 3.0 * y;
            grad[2] = 1.0 - z;
            break;
          case 6:
            grad[0] = 1.0;
            grad[1] = 0;
            grad[2] = 0;
            break;
          case 7:
            grad[0] = 6.0 * z - 1.0;
            grad[1] = 0;
            grad[2] = 6.0 * x;
            break;
          case 8:
            grad[0] = -y * ratio;
            grad[1] = -x * ratio;
            grad[2] = x * y * std::pow(ratio, 2);
            break;
          case 9:
            grad[0] = y * (1.0 - 6.0 * z) * ratio;
            grad[1] = x * (1.0 - 6.0 * z) * ratio;
            grad[2] = 5 * x * y * std::pow(ratio, 2);
            break;
          case 10:
            grad[0] = (-3.0 * std::pow(y, 2) + std::pow(z, 2) - 2.0 * z + 1) *
                      0.5 * ratio;
            grad[1] = -3.0 * x * y * ratio;
            grad[2] = x *
                      (3.0 * std::pow(y, 2) + std::pow(z, 2) - 2.0 * z + 1) *
                      0.5 * std::pow(ratio, 2);
            break;
          case 11:
            grad[0] = 3.0 * x;
            grad[1] = 0.0;
            grad[2] = 1.0 - z;
            break;
          case 12:
            grad[0] = -3.0 * x * y * ratio;
            grad[1] = (-3.0 * std::pow(x, 2) + std::pow(z, 2) - 2.0 * z + 1) *
                      0.5 * ratio;
            grad[2] = y *
                      (3.0 * std::pow(x, 2) + std::pow(z, 2) - 2.0 * z + 1.0) *
                      0.5 * std::pow(ratio, 2);
            break;
          case 13:
            {
              double ratio_13 =
                (ratio == 0.0) ?
                  0.0 :
                  std::pow(z, 3) - 3.0 * std::pow(z, 2) + 3.0 * z - 1.0;
              grad[0] =
                3.0 * x *
                (3.0 * std::pow(y, 2) - std::pow(z, 2) + 2.0 * z - 1.0) * 0.5 *
                std::pow(ratio, 2);
              grad[1] =
                3.0 * y *
                (3.0 * std::pow(x, 2) - std::pow(z, 2) + 2.0 * z - 1.0) * 0.5 *
                std::pow(ratio, 2);
              grad[2] =
                (-9.0 * std::pow(x, 2) * std::pow(y, 2) + std::pow(z, 4) -
                 4.0 * std::pow(z, 3) + 6.0 * std::pow(z, 2) - 4.0 * z + 1) *
                0.5 * ratio_13;
            }
            break;
          default:
            DEAL_II_NOT_IMPLEMENTED();
            break;
        }
      for (unsigned int d = 0; d < dim; ++d)
        if (std::fabs(grad[d]) < 1e-14)
          grad[d] = 0.0;

      return grad;
    }
  DEAL_II_NOT_IMPLEMENTED();
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

  std::cout << "VDM" << std::endl;
  for (unsigned i = 0; i < support_points.size(); ++i)
    {
      for (unsigned j = 0; j < support_points.size(); ++j)
        std::cout << VDM[i][j] << " ";
      std::cout << std::endl;
    }


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
  // AssertIndexRange(this->degree(), 3);
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
  // AssertIndexRange(this->degree(), 3);

  Tensor<1, dim> grad;

  if (false) //(this->degree() == 1)
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
  // else // if (this->degree() == 2)

  {
    AssertIndexRange(i, VDM_inv.m());

    for (unsigned int j = 0; j < VDM_inv.n(); ++j)
      grad += VDM_inv[i][j] * this->compute_jacobi_deriv(j, p);

    for (unsigned int d = 0; d < dim; ++d)
      if (std::fabs(grad[d]) < 1e-14)
        grad[d] = 0.0;
  }

  // else
  //   DEAL_II_NOT_IMPLEMENTED();



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
