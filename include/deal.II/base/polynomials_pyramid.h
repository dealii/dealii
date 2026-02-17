// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2025 by the deal.II authors
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

#include <deal.II/base/point.h>
#include <deal.II/base/scalar_polynomials_base.h>
#include <deal.II/base/tensor.h>

#include <deal.II/lac/full_matrix.h>

DEAL_II_NAMESPACE_OPEN

/**
 * Polynomials defined on pyramid entities. This class is basis of
 * FE_PyramidP.
 * The polynomials are based on @cite Bergot2010. We first use the
 * Jacobi polynomials to construct a modal basis (Proposition 1.10). With the
 * modal basis a Vandermonde matrix is calculated which leads to a nodal basis.
 * For computing the values of the nodal basis the Vandermonde matrix is
 * multiplied with the modal basis vector evaluated at the evaluation point.
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
   * This constructor only works for linear elements.
   */
  ScalarLagrangePolynomialPyramid(const unsigned int degree);

  /*
   * Constructor taking the polynomial @p degree, the number of polynomials
   * @p n_dofs and the support points as input.
   */
  ScalarLagrangePolynomialPyramid(
    const unsigned int             degree,
    const unsigned int             n_dofs,
    const std::vector<Point<dim>> &support_points);

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

  /**
   * @copydoc ScalarPolynomialsBase::compute_value()
   */
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
  /**
   * The Vandermonde matrix evaluates each modal basis function at the chosen
   * nodal points.
   * Applying the inverse of the Vandermonde matrix transforms from the modal
   * basis to the nodal basis.
   */
  FullMatrix<double> vandermonde_matrix_inverse;

  /**
   * Evaluate the orthogonal basis at point @p p. The indices @p i, @p j
   * and @p k corresponde to the polynomial degrees of the Jacobi polynomials,
   * see @cite Bergot2010 proposition 1.10.
   */
  double
  evaluate_orthogonal_basis_function_by_degree(const unsigned int i,
                                               const unsigned int j,
                                               const unsigned int k,
                                               const Point<dim>  &p) const;

  /**
   * Evaluate the orthogonal basis function @p i at point @p p.
   * This function determines the corresponding indices for the Jacobi
   * polynomials and calls the function taking all indices as arguments.
   */
  double
  evaluate_orthogonal_basis_function(const unsigned int i,
                                     const Point<dim>  &p) const;

  /**
   * Evaluate the derivative of the orthogonal basis at point @p p.
   * The indices @p i, @p j and @p k corresponde to the polynomial degrees of
   * the Jacobi polynomials, see @cite Bergot2010 proposition 1.10.
   */
  Tensor<1, dim>
  evaluate_orthogonal_basis_derivative_by_degree(const unsigned int i,
                                                 const unsigned int j,
                                                 const unsigned int k,
                                                 const Point<dim>  &p) const;

  /**
   * Evaluate the derivative of the orthogonal basis function @p i at point
   * @p p. This function determines the corresponding indices for the Jacobi
   * polynomials and calls the function taking all indices as arguments.
   */
  Tensor<1, dim>
  evaluate_orthogonal_basis_derivative(const unsigned int i,
                                       const Point<dim>  &p) const;
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
