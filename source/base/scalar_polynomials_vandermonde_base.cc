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
#include <deal.II/base/scalar_polynomials_vandermonde_base.h>

#include <deal.II/lac/householder.h>
#include <deal.II/lac/vector.h>


DEAL_II_NAMESPACE_OPEN

template <int dim>
ScalarPolynomialsVandermondeBase<dim>::ScalarPolynomialsVandermondeBase(
  const unsigned int degree,
  const unsigned int n_dofs)
  : ScalarPolynomialsBase<dim>(degree, n_dofs)
{}



template <int dim>
void
ScalarPolynomialsVandermondeBase<dim>::reinit(
  const std::vector<Point<dim>> &support_points)
{
  AssertThrow(
    support_points.size() == this->n(),
    ExcNotImplemented(
      "The number of DoFs must be equal to the number of support points."));

  // fill VDM Matrix
  const unsigned int n = support_points.size();
  FullMatrix<double> VDM(n);

  for (unsigned i = 0; i < n; ++i)
    for (unsigned j = 0; j < n; ++j)
      VDM[i][j] = evaluate_orthogonal_basis_function(i, support_points[j]);

  // get the inverse matrix
  Householder<double> householder(VDM);
  Vector<double>      e(n);
  Vector<double>      x(n);

  // assume e_j are unit vectors
  // compute (VDM^-1)_ij by using (VDM^-1)_ij = e_i (VDM^-1) e_j
  // loop over all vectors e_j
  for (unsigned int j = 0; j < n; ++j)
    {
      e    = 0.;
      e[j] = 1.;

      x = 0.;
      // get x = (VDM^-1) e_j
      householder.least_squares(x, e);

      for (unsigned int i = 0; i < n; ++i)
        // (VDM^-1)_ij = e_i (VDM^-1) e_j
        VDM(i, j) = x[i];
    }

  vandermonde_matrix_inverse = VDM;

  // clean up small values
  for (unsigned int i = 0; i < n; ++i)
    for (unsigned int j = 0; j < n; ++j)
      if (std::abs(vandermonde_matrix_inverse[i][j]) < 1e-14)
        vandermonde_matrix_inverse[i][j] = 0.;
}



template <int dim>
double
ScalarPolynomialsVandermondeBase<dim>::compute_value(const unsigned int i,
                                                     const Point<dim>  &p) const
{
  AssertIndexRange(i, vandermonde_matrix_inverse.m());

  double result = 0.;
  for (unsigned int j = 0; j < vandermonde_matrix_inverse.n(); ++j)
    result += vandermonde_matrix_inverse[i][j] *
              evaluate_orthogonal_basis_function(j, p);

  if (std::fabs(result) < 1e-14)
    result = 0.0;

  return result;
}



template <int dim>
Tensor<1, dim>
ScalarPolynomialsVandermondeBase<dim>::compute_grad(const unsigned int i,
                                                    const Point<dim>  &p) const
{
  AssertIndexRange(i, vandermonde_matrix_inverse.m());

  Tensor<1, dim> grad;

  grad = 0.;
  for (unsigned int j = 0; j < vandermonde_matrix_inverse.n(); ++j)
    grad += vandermonde_matrix_inverse[i][j] *
            evaluate_orthogonal_basis_derivative(j, p);

  if constexpr (dim > 0)
    for (unsigned int d = 0; d < dim; ++d)
      if (std::fabs(grad[d]) < 1e-14)
        grad[d] = 0.0;

  return grad;
}



template <int dim>
Tensor<2, dim>
ScalarPolynomialsVandermondeBase<dim>::compute_grad_grad(
  const unsigned int i,
  const Point<dim>  &p) const
{
  AssertIndexRange(i, vandermonde_matrix_inverse.m());

  Tensor<2, dim> grad_grad;

  grad_grad = 0.;
  for (unsigned int j = 0; j < vandermonde_matrix_inverse.n(); ++j)
    grad_grad += vandermonde_matrix_inverse[i][j] *
                 evaluate_orthogonal_basis_2nd_derivative(j, p);

  if constexpr (dim > 0)
    for (unsigned int d = 0; d < dim; ++d)
      for (unsigned int e = 0; e < dim; ++e)
        if (std::fabs(grad_grad[d][e]) < 1e-14)
          grad_grad[d][e] = 0.0;

  return grad_grad;
}



template <int dim>
void
ScalarPolynomialsVandermondeBase<dim>::evaluate(
  const Point<dim>            &unit_point,
  std::vector<double>         &values,
  std::vector<Tensor<1, dim>> &grads,
  std::vector<Tensor<2, dim>> &grad_grads,
  std::vector<Tensor<3, dim>> &third_derivatives,
  std::vector<Tensor<4, dim>> &fourth_derivatives) const
{
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
ScalarPolynomialsVandermondeBase<dim>::compute_1st_derivative(
  const unsigned int i,
  const Point<dim>  &p) const
{
  return compute_grad(i, p);
}



template <int dim>
Tensor<2, dim>
ScalarPolynomialsVandermondeBase<dim>::compute_2nd_derivative(
  const unsigned int i,
  const Point<dim>  &p) const
{
  return compute_grad_grad(i, p);
}



template <int dim>
Tensor<3, dim>
ScalarPolynomialsVandermondeBase<dim>::compute_3rd_derivative(
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
ScalarPolynomialsVandermondeBase<dim>::compute_4th_derivative(
  const unsigned int i,
  const Point<dim>  &p) const
{
  (void)i;
  (void)p;

  DEAL_II_NOT_IMPLEMENTED();

  return {};
}



template <int dim>
Tensor<2, dim>
ScalarPolynomialsVandermondeBase<dim>::
  evaluate_orthogonal_basis_2nd_derivative_by_degree(const unsigned int i,
                                                     const unsigned int j,
                                                     const unsigned int k,
                                                     const Point<dim>  &p) const
{
  (void)i;
  (void)j;
  (void)k;
  (void)p;

  DEAL_II_NOT_IMPLEMENTED();

  return {};
}



template <int dim>
Tensor<2, dim>
ScalarPolynomialsVandermondeBase<dim>::evaluate_orthogonal_basis_2nd_derivative(
  const unsigned int i,
  const Point<dim>  &p) const
{
  (void)i;
  (void)p;

  DEAL_II_NOT_IMPLEMENTED();

  return {};
}



template class ScalarPolynomialsVandermondeBase<0>;
template class ScalarPolynomialsVandermondeBase<1>;
template class ScalarPolynomialsVandermondeBase<2>;
template class ScalarPolynomialsVandermondeBase<3>;

DEAL_II_NAMESPACE_CLOSE
