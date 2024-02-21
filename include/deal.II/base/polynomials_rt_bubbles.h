// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_polynomials_rt_bubbles_h
#define dealii_polynomials_rt_bubbles_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/point.h>
#include <deal.II/base/polynomial.h>
#include <deal.II/base/polynomials_raviart_thomas.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/tensor_polynomials_base.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN

/**
 * This class implements the <i>H<sup>div</sup></i>-conforming, vector-valued
 * enhanced Raviart-Thomas polynomials.
 *
 * Similarly to the classical Raviart-Thomas space, the enhanced Raviart-Thomas
 * polynomials are constructed such that the divergence is in the tensor product
 * polynomial space <i>Q<sub>k-1</sub></i>.
 *
 * This space is of the form <i>V<sub>k</sub> = RT<sub>k-1</sub> +
 * B<sub>k</sub></i>, where <i>B<sub>k</sub></i> is defined as follows:
 * <dl>
 * <dt> In 2d:</dt>
 * <dd>
 * @f{align*}{
 *  B_k^1(E) = \text{span}\left\{x^{a_1-1} y^{a_2}\begin{pmatrix} (a_2+1) x \\
 *    -a_1 y \end{pmatrix}\text{ : } a_2=k \right\} \\
 *  B_k^2(E) = \text{span}\left\{x^{b_1} y^{b_2-1}\begin{pmatrix} -b_2 x \\
 *     (b_1+1) y \end{pmatrix}\text{ : } b_1=k \right\}
 * @f}
 * </dd>
 *
 * <dt> In 3d: </dt>
 * <dd>
 *  @f{align*}{
 *   B_k^1(E) = \text{span}\left\{x^{a_1-1} y^{a_2} z^{a_3}\begin{pmatrix}
 * (a_2+a_3+2) x \\
 *     -a_1 y \\ -a_1 z \end{pmatrix}\text{ : } a_2=k \text{ or } a_3=k
 * \right\},\\
 *   B_k^2(E) = \text{span}\left\{x^{b_1} y^{b_2-1} z^{b_3}\begin{pmatrix} -b_2
 * x \\
 *     (b_1+b_3+2) y \\ -b_2 z \end{pmatrix}\text{ : } b_1=k \text{ or } b_3=k
 * \right\},\\
 *   B_k^3(E) = \text{span}\left\{x^{c_1}y^{c_2}z^{c_3-1}\begin{pmatrix} -c_3 x
 * \\ -c_3y \\ (c_1+c_2+2)z \end{pmatrix}\text{ : } c_1=k \text{ or } c_2=k
 * \right\},
 *  @f}
 * </dd>
 * </dl>
 * where $0 \le a_1, a_2, a_3 \le k$.
 *
 * @note Unlike the classical Raviart-Thomas space, the lowest order for the enhanced space
 * is 1, similarly to the Brezzi-Douglas-Marini (BDM) polynomial space.
 *
 * The total dimension of the space <i>dim(V<sub>k</sub>) = d*(k+1)^d</i>, where
 * <i>d</i> is the space dimension. This allows to associate shape functions
 * with the Gauss-Lobatto quadrature points as shown in the figures below.
 *
 * <table> <tr> <td align="center">
 * @image html rtbubbles.png
 * </td></tr>
 *
 * <tr> <td align="center"> Left - $2d,\,k=3$,
 * right - $3d,\,k=2$.</td></tr> </table>
 *
 * @ingroup Polynomials
 */

template <int dim>
class PolynomialsRT_Bubbles : public TensorPolynomialsBase<dim>
{
public:
  /**
   * Constructor. Creates all basis functions for RT_bubbles polynomials of
   * given degree.
   */
  PolynomialsRT_Bubbles(const unsigned int k);

  /**
   * Computes the value and the first and second derivatives of each RT_bubbles
   * polynomial at @p unit_point.
   *
   * The size of the vectors must either be zero or equal <tt>n()</tt>.  In
   * the first case, the function will not compute these values.
   */
  void
  evaluate(const Point<dim>            &unit_point,
           std::vector<Tensor<1, dim>> &values,
           std::vector<Tensor<2, dim>> &grads,
           std::vector<Tensor<3, dim>> &grad_grads,
           std::vector<Tensor<4, dim>> &third_derivatives,
           std::vector<Tensor<5, dim>> &fourth_derivatives) const override;

  /**
   * Return the name of the space, which is <tt>RT_Bubbles</tt>.
   */
  std::string
  name() const override;

  /**
   * Return the number of polynomials in the space <tt>RT_Bubbles(degree)</tt>
   * without requiring to build an object of PolynomialsRT-Bubbles. This is
   * required by the FiniteElement classes.
   */
  static unsigned int
  n_polynomials(const unsigned int degree);

  /**
   * @copydoc TensorPolynomialsBase::clone()
   */
  virtual std::unique_ptr<TensorPolynomialsBase<dim>>
  clone() const override;

private:
  /**
   * An object representing the Raviart-Thomas part of the space
   */
  const PolynomialsRaviartThomas<dim> raviart_thomas_space;

  /**
   * Storage for monomials, we need all polynomials from degree zero
   * to <i>k+1</i>.
   */
  std::vector<Polynomials::Polynomial<double>> monomials;
};


template <int dim>
inline std::string
PolynomialsRT_Bubbles<dim>::name() const
{
  return "RT_bubbles";
}


DEAL_II_NAMESPACE_CLOSE

#endif
