// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#ifndef dealii_polynomials_nedelec_nodal_h
#define dealii_polynomials_nedelec_nodal_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/point.h>
#include <deal.II/base/polynomial.h>
#include <deal.II/base/polynomial_space.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/tensor_polynomials_base.h>
#include <deal.II/base/tensor_product_polynomials.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN

/**
 * This class implements the first family <i>H<sup>curl</sup></i>-conforming,
 * vector-valued polynomials, proposed by J.-C. Nédélec in 1980 (Numer.
 * Math. 35).
 *
 * The Nédélec polynomials are constructed such that the curl is in the
 * tensor product polynomial space <i>Q<sub>k</sub></i>. Therefore, the
 * polynomial order of each component must be one order higher in the
 * corresponding two directions, yielding the polynomial spaces
 * <i>(Q<sub>k,k+1</sub>, Q<sub>k+1,k</sub>)</i> and
 * <i>(Q<sub>k,k+1,k+1</sub>, Q<sub>k+1,k,k+1</sub>,
 * Q<sub>k+1,k+1,k</sub>)</i> in 2d and 3d, resp.
 *
 * @ingroup Polynomials
 */
template <int dim>
class PolynomialsNedelecNodal : public TensorPolynomialsBase<dim>
{
public:
  /**
   * Constructor. Creates all basis functions for Nédélec polynomials of
   * given degree.
   *
   * @arg k: the degree of the Nédélec space, which is the degree of the
   * largest tensor product polynomial space <i>Q<sub>k</sub></i> contained.
   */
  PolynomialsNedelecNodal(const unsigned int degree);


  /**
   * Compute the value and the first and second derivatives of each Nedelec
   * polynomial at @p unit_point.
   *
   * The size of the vectors must either be zero or equal <tt>n()</tt>.  In
   * the first case, the function will not compute these values.
   *
   * If you need values or derivatives of all tensor product polynomials then
   * use this function, rather than using any of the <tt>compute_value</tt>,
   * <tt>compute_grad</tt> or <tt>compute_grad_grad</tt> functions, see below,
   * in a loop over all tensor product polynomials.
   */
  void
  evaluate(const Point<dim>            &unit_point,
           std::vector<Tensor<1, dim>> &values,
           std::vector<Tensor<2, dim>> &grads,
           std::vector<Tensor<3, dim>> &grad_grads,
           std::vector<Tensor<4, dim>> &third_derivatives,
           std::vector<Tensor<5, dim>> &fourth_derivatives) const override;

  /**
   * Return the name of the space, which is <tt>Nedelec</tt>.
   */
  std::string
  name() const override;

  /**
   * Return the number of polynomials in the space <tt>N(degree)</tt> without
   * requiring to build an object of PolynomialsNedelec. This is required by
   * the FiniteElement classes.
   */
  static unsigned int
  n_polynomials(const unsigned int degree);

  /**
   * Compute the lexicographic to hierarchic numbering underlying this class,
   * computed as a free function.
   */
  static std::vector<unsigned int>
  get_lexicographic_numbering(const unsigned int degree);

  /**
   * @copydoc TensorPolynomialsBase::clone()
   */
  virtual std::unique_ptr<TensorPolynomialsBase<dim>>
  clone() const override;

  /**
   * Compute the generalized support points in the ordering used by the
   * polynomial shape functions. Note that these points are not support points
   * in the classical sense as the Lagrange polynomials of the different
   * components have different points, which need to be combined in terms of
   * Piola transforms.
   */
  std::vector<Point<dim>>
  get_polynomial_support_points() const;

private:
  /**
   * An object representing the polynomial space for a single component. We
   * can re-use it by rotating the coordinates of the evaluation point.
   */
  const AnisotropicPolynomials<dim> polynomial_space;

  /**
   * The given degree
   */
  const unsigned int deg;

  /**
   * Renumbering from lexicographic to hierarchic order.
   */
  std::vector<unsigned int> renumber_lexicographic_to_hierarchic;

  /**
   * Renumbering from hierarchic to lexicographic order. Inverse of
   * renumber_lexicographic_to_hierarchic.
   */
  std::vector<unsigned int> renumber_hierarchic_to_lexicographic;

  /**
   * Renumbering from shifted polynomial spaces to lexicographic one
   */
  std::array<std::vector<unsigned int>, dim> renumber_aniso;

  /**
   * A static member function that creates the polynomial space we use to
   * initialize the #polynomial_space member variable.
   */
  static std::vector<std::vector<Polynomials::Polynomial<double>>>
  create_polynomials(const unsigned int degree);
};

DEAL_II_NAMESPACE_CLOSE

#endif
