// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#ifndef dealii_polynomials_pyramid_h
#define dealii_polynomials_pyramid_h

#include <deal.II/base/config.h>

#include <deal.II/base/scalar_polynomials_base.h>

#include <deal.II/lac/full_matrix.h>

DEAL_II_NAMESPACE_OPEN

/**
 * Polynomials defined on pyramid entities. This class is basis of
 * FE_PyramidP.
 */
template <int dim>
class ScalarLagrangePolynomialPyramid : public ScalarPolynomialsBase<dim>
{
public:
  /**
   * Make the dimension available to the outside.
   */
  static constexpr unsigned int dimension = dim;

  /*
   * Constructor taking the polynomial @p degree as input.
   */
  ScalarLagrangePolynomialPyramid(const unsigned int            degree,
                                  const unsigned int            n_dofs,
                                  const std::vector<Point<dim>> support_points);

  /**
   * @copydoc ScalarPolynomialsBase::evaluate()
   *
   * @note Currently, only the vectors @p values and @p grads are filled.
   */
  void
  evaluate(const Point<dim>            &unit_point,
           std::vector<double>         &values,
           std::vector<Tensor<1, dim>> &grads,
           std::vector<Tensor<2, dim>> &grad_grads,
           std::vector<Tensor<3, dim>> &third_derivatives,
           std::vector<Tensor<4, dim>> &fourth_derivatives) const override;

  double
  compute_value(const unsigned int i, const Point<dim> &p) const override;

  /**
   * @copydoc ScalarPolynomialsBase::compute_derivative()
   *
   * @note Currently, only implemented for first derivative.
   */
  template <int order>
  Tensor<order, dim>
  compute_derivative(const unsigned int i, const Point<dim> &p) const;

  Tensor<1, dim>
  compute_1st_derivative(const unsigned int i,
                         const Point<dim>  &p) const override;

  Tensor<2, dim>
  compute_2nd_derivative(const unsigned int i,
                         const Point<dim>  &p) const override;

  /**
   * @copydoc ScalarPolynomialsBase::compute_3rd_derivative()
   *
   * @note Not implemented yet.
   */
  Tensor<3, dim>
  compute_3rd_derivative(const unsigned int i,
                         const Point<dim>  &p) const override;

  /**
   * @copydoc ScalarPolynomialsBase::compute_4th_derivative()
   *
   * @note Not implemented yet.
   */
  Tensor<4, dim>
  compute_4th_derivative(const unsigned int i,
                         const Point<dim>  &p) const override;

  /**
   * @copydoc ScalarPolynomialsBase::compute_grad()
   *
   * @note Not implemented yet.
   */
  Tensor<1, dim>
  compute_grad(const unsigned int i, const Point<dim> &p) const override;

  /**
   * @copydoc ScalarPolynomialsBase::compute_grad_grad()
   *
   * @note Not implemented yet.
   */
  Tensor<2, dim>
  compute_grad_grad(const unsigned int i, const Point<dim> &p) const override;

  std::string
  name() const override;

  virtual std::unique_ptr<ScalarPolynomialsBase<dim>>
  clone() const override;

private:
  FullMatrix<double> VDM_inv;

  double
  compute_jacobi_basis_functions(const unsigned int i,
                                 const Point<dim>  &p) const;

  Tensor<1, dim>
  compute_jacobi_deriv_basis_functions(const unsigned int i,
                                       const Point<dim>  &p) const;

  double
  compute_polynomial_space(const unsigned int i,
                           const unsigned int j,
                           const unsigned int k,
                           const Point<dim>  &p) const;

  double
  compute_jacobi_basis(const unsigned int i, const Point<dim> &p) const;

  Tensor<1, dim>
  compute_polynomial_space_deriv(const unsigned int i,
                                 const unsigned int j,
                                 const unsigned int k,
                                 const Point<dim>  &p) const;

  Tensor<1, dim>
  compute_jacobi_deriv(const unsigned int i, const Point<dim> &p) const;
};



template <int dim>
template <int order>
Tensor<order, dim>
ScalarLagrangePolynomialPyramid<dim>::compute_derivative(
  const unsigned int i,
  const Point<dim>  &p) const
{
  Tensor<order, dim> der;

  Assert(order == 1, ExcNotImplemented());
  const auto grad = compute_grad(i, p);

  for (unsigned int i = 0; i < dim; ++i)
    der[i] = grad[i];

  return der;
}

DEAL_II_NAMESPACE_CLOSE

#endif
