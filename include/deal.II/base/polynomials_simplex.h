// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2021 - 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


#ifndef dealii_polynomials_simplex_h
#define dealii_polynomials_simplex_h

#include <deal.II/base/config.h>

#include <deal.II/base/point.h>
#include <deal.II/base/scalar_polynomials_vandermonde_base.h>
#include <deal.II/base/tensor.h>

#include <deal.II/lac/full_matrix.h>

DEAL_II_NAMESPACE_OPEN

/**
 * Polynomials defined on simplex entities. This class can be a basis of
 * FE_SimplexP.
 * We first use the Jacobi polynomials given in @cite Hesthaven2007 to construct
 * a modal basis. With the modal basis a Vandermonde matrix is calculated which
 * leads to a nodal basis. For computing the values of the nodal basis the
 * Vandermonde matrix is multiplied with the modal basis vector evaluated at the
 * evaluation point.
 */
template <int dim>
class ScalarLagrangePolynomialSimplex
  : public ScalarPolynomialsVandermondeBase<dim>
{
public:
  /**
   * Make the dimension available to the outside.
   */
  static constexpr unsigned int dimension = dim;

  /*
   * Constructor taking the polynomial @p degree, the number of polynomials
   * @p n_dofs and the support points as input.
   */
  ScalarLagrangePolynomialSimplex(
    const unsigned int             degree,
    const std::vector<Point<dim>> &support_points);

  std::string
  name() const override;

  virtual std::unique_ptr<ScalarPolynomialsBase<dim>>
  clone() const override;

private:
  /**
   * Evaluate the orthogonal basis at point @p p. The indices @p i, @p j
   * and @p k correspond to the polynomial degrees of the Jacobi polynomials
   * given in @cite Hesthaven2007.
   */
  double
  evaluate_orthogonal_basis_function_by_degree(
    const unsigned int i,
    const unsigned int j,
    const unsigned int k,
    const Point<dim>  &p) const override;

  /**
   * Evaluate the orthogonal basis function @p i at point @p p.
   * This function determines the corresponding indices for the Jacobi
   * polynomials and calls the function taking all indices as arguments.
   */
  double
  evaluate_orthogonal_basis_function(const unsigned int i,
                                     const Point<dim>  &p) const override;

  /**
   * Evaluate the derivative of the orthogonal basis at point @p p.
   * The indices @p i, @p j and @p k correspond to the polynomial degrees of
   * the Jacobi polynomials given in  @cite Hesthaven2007.
   */
  Tensor<1, dim>
  evaluate_orthogonal_basis_derivative_by_degree(
    const unsigned int i,
    const unsigned int j,
    const unsigned int k,
    const Point<dim>  &p) const override;

  /**
   * Evaluate the derivative of the orthogonal basis function @p i at point
   * @p p. This function determines the corresponding indices for the Jacobi
   * polynomials and calls the function taking all indices as arguments.
   */
  Tensor<1, dim>
  evaluate_orthogonal_basis_derivative(const unsigned int i,
                                       const Point<dim>  &p) const override;
};


DEAL_II_NAMESPACE_CLOSE

#endif
