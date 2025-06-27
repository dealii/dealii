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


#ifndef dealii_base_polynomials_wedge_h
#define dealii_base_polynomials_wedge_h

#include <deal.II/base/config.h>

#include <deal.II/base/ndarray.h>
#include <deal.II/base/polynomials_barycentric.h>
#include <deal.II/base/scalar_polynomials_base.h>

DEAL_II_NAMESPACE_OPEN


namespace internal
{
  /**
   * Decompose the shape-function index of a linear wedge into an index
   * to access the right shape function within the triangle and within
   * the line.
   */
  constexpr dealii::ndarray<unsigned int, 6, 2> wedge_table_1{
    {{{0, 0}}, {{1, 0}}, {{2, 0}}, {{0, 1}}, {{1, 1}}, {{2, 1}}}};

  /**
   * Decompose the shape-function index of a quadratic wedge into an index
   * to access the right shape function within the triangle and within
   * the line.
   */
  constexpr dealii::ndarray<unsigned int, 18, 2> wedge_table_2{{{{0, 0}},
                                                                {{1, 0}},
                                                                {{2, 0}},
                                                                {{0, 1}},
                                                                {{1, 1}},
                                                                {{2, 1}},
                                                                {{3, 0}},
                                                                {{4, 0}},
                                                                {{5, 0}},
                                                                {{3, 1}},
                                                                {{4, 1}},
                                                                {{5, 1}},
                                                                {{0, 2}},
                                                                {{1, 2}},
                                                                {{2, 2}},
                                                                {{3, 2}},
                                                                {{4, 2}},
                                                                {{5, 2}}}};
} // namespace internal


/**
 * Polynomials defined on wedge entities. This class is basis of
 * FE_WedgeP.
 *
 * The polynomials are created via a tensor product of a
 * BarycentricPolynomials<2>::get_fe_p_basis(degree) and a
 * BarycentricPolynomials<1>::get_fe_p_basis(degree), however, are
 * re-numerated to better match the definition of FiniteElement.
 */
template <int dim>
class ScalarLagrangePolynomialWedge : public ScalarPolynomialsBase<dim>
{
public:
  /**
   * Make the dimension available to the outside.
   */
  static constexpr unsigned int dimension = dim;

  /*
   * Constructor taking the polynomial @p degree as input.
   *
   * @note Currently, only linear (degree=1) and quadratic polynomials
   *   (degree=2) are implemented.
   */
  ScalarLagrangePolynomialWedge(const unsigned int degree);

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
   * Scalar polynomials defined on a triangle.
   */
  const BarycentricPolynomials<2> poly_tri;

  /**
   * Scalar polynomials defined on a line.
   */
  const BarycentricPolynomials<1> poly_line;
};



template <int dim>
template <int order>
Tensor<order, dim>
ScalarLagrangePolynomialWedge<dim>::compute_derivative(
  const unsigned int i,
  const Point<dim>  &p) const
{
  Tensor<order, dim> der;

  AssertDimension(order, 1);
  const auto grad = compute_grad(i, p);

  for (unsigned int i = 0; i < dim; ++i)
    der[i] = grad[i];

  return der;
}

DEAL_II_NAMESPACE_CLOSE

#endif
