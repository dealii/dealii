// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2004 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_polynomials_BDM_h
#define dealii_polynomials_BDM_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/mutex.h>
#include <deal.II/base/point.h>
#include <deal.II/base/polynomial.h>
#include <deal.II/base/polynomial_space.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/tensor_polynomials_base.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN

/**
 * This class implements the <i>H<sup>div</sup></i>-conforming, vector-valued
 * Brezzi-Douglas-Marini (<i> BDM </i>) polynomials described in Brezzi and
 * Fortin's <i>Mixed and Hybrid Finite Element Methods</i> (refer to pages
 * 119 - 124).
 *
 * The <i> BDM </i> polynomial space contain the entire  $(P_{k})^{n}$
 * space (constructed with PolynomialSpace Legendre polynomials) as well as
 * part of $(P_{k+1})^{n}$
 * (ie. $(P_{k})^{n} \subset BDM_{k} \subset (P_{k+1})^{n}$).  Furthermore,
 * $BDM_{k}$ elements are designed so that
 * $\nabla \cdot q \in P_{k-1} (K)$ and $q \cdot n |_{e_{i}} \in P_{k}(e_{i})$.
 * More details
 * of two and three dimensional $BDM_{k}$ elements are given below.
 *<dl>
 *   <dt> In 2d:
 *   <dd> $ BDM_{k} = \{\mathbf{q} | \mathbf{q} = p_{k} (x,y) +
 *      r \; \text{curl} (x^{k+1}y) + s \;
 *      \text{curl} (xy^{k+1}), p_{k} \in (P_{k})^{2} \}$.
 *
 *   Note: the curl of a scalar function is given by $\text{curl}(f(x,y)) =
 *   \begin{pmatrix} f_{y}(x,y) \\ -f_{x}(x,y) \end{pmatrix}$.
 *
 *   The basis used to construct the $BDM_{1}$ shape functions is
 *   @f{align*}{
 *     \phi_0 = \begin{pmatrix} 1 \\ 0 \end{pmatrix},
 *     \phi_1 = \begin{pmatrix} -\sqrt{3}+2\sqrt{3}x \\ 0 \end{pmatrix},
 *     \phi_2 = \begin{pmatrix} -\sqrt{3}+2\sqrt{3}y \\ 0 \end{pmatrix},
 *     \phi_3 = \begin{pmatrix} 0 \\ 1 \end{pmatrix},
 *     \phi_4 = \begin{pmatrix} 0 \\ -\sqrt{3}+2\sqrt{3}x \end{pmatrix},
 *     \phi_5 = \begin{pmatrix} 0 \\ -\sqrt{3}+2\sqrt{3}y \end{pmatrix},
 *     \phi_6 = \begin{pmatrix} x^2 \\ -2xy \end{pmatrix},
 *     \phi_7 = \begin{pmatrix} 2xy \\ -y^2 \end{pmatrix}.
 *   @f}
 *
 *   The dimension of the $BDM_{k}$ space is
 * $(k+1)(k+2)+2$, with $k+1$ unknowns per
 * edge and $k(k-1)$ interior unknowns.
 *
 *   <dt> In 3d:
 *   <dd> $ BDM_{k} =
 *        \{\mathbf{q} | \mathbf{q} = p_{k} (x,y,z)
 *        + \sum_{i=0}^{k} (
 *        r_{i} \; \text{curl}
 *        \begin{pmatrix} 0\\0\\xy^{i+1}z^{k-i} \end{pmatrix}
 *        + s_{i} \; \text{curl}
 *        \begin{pmatrix} yz^{i+1}x^{k-i}\\0\\0 \end{pmatrix}
 *        + t_{i} \; \text{curl}
 *        \begin{pmatrix}0\\zx^{i+1}y^{k-i}\\0\end{pmatrix})
 *        , p_{k} \in (P_{k})^{3} \}$.
 *
 *   Note: the 3d description of $BDM_{k}$ is not unique.  See <i>Mixed and
 *   Hybrid Finite Element Methods</i> page 122 for an alternative definition.
 *
 *   The dimension of the $BDM_{k}$ space is
 *   $\dfrac{(k+1)(k+2)(k+3)}{2}+3(k+1)$, with $\dfrac{(k+1)(k+2)}{2}$
 *   unknowns per face and $\dfrac{(k-1)k(k+1)}{2}$ interior unknowns.
 *
 *</dl>
 *
 *
 *
 * @ingroup Polynomials
 */
template <int dim>
class PolynomialsBDM : public TensorPolynomialsBase<dim>
{
public:
  /**
   * Constructor. Creates all basis functions for BDM polynomials of given
   * degree.
   *
   * @arg k: the degree of the BDM-space, which is the degree of the largest
   * complete polynomial space <i>P<sub>k</sub></i> contained in the BDM-
   * space.
   */
  PolynomialsBDM(const unsigned int k);

  /**
   * Compute the value and the first and second derivatives of each BDM
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
   * Return the name of the space, which is <tt>BDM</tt>.
   */
  std::string
  name() const override;

  /**
   * Return the number of polynomials in the space <tt>BDM(degree)</tt>
   * without requiring to build an object of PolynomialsBDM. This is required
   * by the FiniteElement classes.
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
   * An object representing the polynomial space used here. The constructor
   * fills this with the monomial basis.
   */
  const PolynomialSpace<dim> polynomial_space;

  /**
   * Storage for monomials. In 2d, this is just the polynomial of order
   * <i>k</i>. In 3d, we need all polynomials from degree zero to <i>k</i>.
   */
  std::vector<Polynomials::Polynomial<double>> monomials;

  /**
   * A mutex that guards the following scratch arrays.
   */
  mutable Threads::Mutex mutex;

  /**
   * Auxiliary memory.
   */
  mutable std::vector<double> p_values;

  /**
   * Auxiliary memory.
   */
  mutable std::vector<Tensor<1, dim>> p_grads;

  /**
   * Auxiliary memory.
   */
  mutable std::vector<Tensor<2, dim>> p_grad_grads;

  /**
   * Auxiliary memory.
   */
  mutable std::vector<Tensor<3, dim>> p_third_derivatives;

  /**
   * Auxiliary memory.
   */
  mutable std::vector<Tensor<4, dim>> p_fourth_derivatives;
};


template <int dim>
inline std::string
PolynomialsBDM<dim>::name() const
{
  return "BDM";
}


DEAL_II_NAMESPACE_CLOSE

#endif
