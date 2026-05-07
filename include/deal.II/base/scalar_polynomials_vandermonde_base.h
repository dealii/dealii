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

#ifndef dealii_scalar_polynomials_vandermonde_base_h
#define dealii_scalar_polynomials_vandermonde_base_h


#include <deal.II/base/config.h>

#include <deal.II/base/point.h>
#include <deal.II/base/scalar_polynomials_base.h>
#include <deal.II/base/tensor.h>

#include <deal.II/lac/full_matrix.h>

DEAL_II_NAMESPACE_OPEN

/**
 * This class provides a framework for finite elements using a nodal polynomial
 * basis, where the polynomial basis is constructed by evaluating the values of
 * a modal basis first that then gets transformed to the actual nodal polynomial
 * values via a Vandermonde matrix. This is a common approach for high-order
 * bases on general point distributions.
 *
 * Any derived class must provide the most basic properties for the modal basis
 * like evaluate_orthogonal_basis_function_by_degree(),
 * evaluate_orthogonal_basis_function(),
 * evaluate_orthogonal_basis_derivative_by_degree() and
 * evaluate_orthogonal_basis_derivative().
 */
template <int dim>
class ScalarPolynomialsVandermondeBase : public ScalarPolynomialsBase<dim>
{
public:
  /**
   * Constructor. This takes the degree @p degree of the space and the number
   * of polynomials @p n.
   */
  ScalarPolynomialsVandermondeBase(const unsigned int degree,
                                   const unsigned int n_dofs);

  /*
   * Virtual destructor. Makes sure that pointers to this class are deleted
   * properly.
   */
  virtual ~ScalarPolynomialsVandermondeBase() = default;

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

  /**
   * @copydoc ScalarPolynomialsBase::compute_1st_derivative()
   */
  Tensor<1, dim>
  compute_1st_derivative(const unsigned int i,
                         const Point<dim>  &p) const override;

  /**
   * @copydoc ScalarPolynomialsBase::compute_2nd_derivative()
   *
   * @note Not implemented yet.
   */
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

protected:
  /**
   * The Vandermonde matrix evaluates each modal basis function at the chosen
   * nodal points.
   * Applying the inverse of the Vandermonde matrix transforms from the modal
   * basis to the nodal basis.
   */
  FullMatrix<double> vandermonde_matrix_inverse;

  /*
   * Evaluate the modal basis at all support points @p support_points to construct
   * the Vandermonde matrix and invert it.
   */
  void
  reinit(const std::vector<Point<dim>> &support_points);

  /**
   * Evaluate the orthogonal basis at point @p p. The indices @p i, @p j
   * and @p k correspond to the polynomial degrees of the Jacobi polynomials.
   */
  virtual double
  evaluate_orthogonal_basis_function_by_degree(const unsigned int i,
                                               const unsigned int j,
                                               const unsigned int k,
                                               const Point<dim>  &p) const = 0;

  /**
   * Evaluate the orthogonal basis function @p i at point @p p.
   * This function determines the corresponding indices for the Jacobi
   * polynomials and calls the function taking all indices as arguments.
   */
  virtual double
  evaluate_orthogonal_basis_function(const unsigned int i,
                                     const Point<dim>  &p) const = 0;

  /**
   * Evaluate the derivative of the orthogonal basis at point @p p.
   * The indices @p i, @p j and @p k correspond to the polynomial degrees of
   * the Jacobi polynomials.
   */
  virtual Tensor<1, dim>
  evaluate_orthogonal_basis_derivative_by_degree(const unsigned int i,
                                                 const unsigned int j,
                                                 const unsigned int k,
                                                 const Point<dim> &p) const = 0;

  /**
   * Evaluate the derivative of the orthogonal basis function @p i at point
   * @p p. This function determines the corresponding indices for the Jacobi
   * polynomials and calls the function taking all indices as arguments.
   */
  virtual Tensor<1, dim>
  evaluate_orthogonal_basis_derivative(const unsigned int i,
                                       const Point<dim>  &p) const = 0;
};



template <int dim>
template <int order>
Tensor<order, dim>
ScalarPolynomialsVandermondeBase<dim>::compute_derivative(
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
