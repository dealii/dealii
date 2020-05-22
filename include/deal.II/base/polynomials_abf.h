// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

#ifndef dealii_polynomials_abf_h
#define dealii_polynomials_abf_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/point.h>
#include <deal.II/base/polynomial.h>
#include <deal.II/base/polynomial_space.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/tensor_polynomials_base.h>
#include <deal.II/base/tensor_product_polynomials.h>
#include <deal.II/base/thread_management.h>

#include <memory>
#include <vector>

DEAL_II_NAMESPACE_OPEN

/**
 * This class implements the <i>H<sup>div</sup></i>-conforming, vector-valued
 * Arnold-Boffi-Falk polynomials as described in the article by Arnold-Boffi-
 * Falk: Quadrilateral H(div) finite elements, SIAM J. Numer. Anal. Vol.42,
 * No.6, pp.2429-2451
 *
 *
 * The ABF polynomials are constructed such that the divergence is in the
 * tensor product polynomial space <i>Q<sub>k</sub></i>. Therefore, the
 * polynomial order of each component must be two orders higher in the
 * corresponding direction, yielding the polynomial spaces
 * <i>(Q<sub>k+2,k</sub>, Q<sub>k,k+2</sub>)</i> and <i>(Q<sub>k+2,k,k</sub>,
 * Q<sub>k,k+2,k</sub>, Q<sub>k,k,k+2</sub>)</i> in 2D and 3D, resp.
 *
 * @ingroup Polynomials
 * @author Oliver Kayser-Herold, based on code from Guido Kanschat
 * @date 2006
 */
template <int dim>
class PolynomialsABF : public TensorPolynomialsBase<dim>
{
public:
  /**
   * Constructor. Creates all basis functions for Raviart-Thomas polynomials
   * of given degree.
   *
   * @arg k: the degree of the Raviart-Thomas-space, which is the degree of
   * the largest tensor product polynomial space <i>Q<sub>k</sub></i>
   * contained.
   */
  PolynomialsABF(const unsigned int k);

  /**
   * Compute the value and the first and second derivatives of each Raviart-
   * Thomas polynomial at @p unit_point.
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
  evaluate(const Point<dim> &           unit_point,
           std::vector<Tensor<1, dim>> &values,
           std::vector<Tensor<2, dim>> &grads,
           std::vector<Tensor<3, dim>> &grad_grads,
           std::vector<Tensor<4, dim>> &third_derivatives,
           std::vector<Tensor<5, dim>> &fourth_derivatives) const override;

  /**
   * Return the name of the space, which is <tt>ABF</tt>.
   */
  std::string
  name() const override;

  /**
   * Return the number of polynomials in the space <tt>RT(degree)</tt> without
   * requiring to build an object of PolynomialsABF. This is required by the
   * FiniteElement classes.
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
   * An object representing the polynomial space for a single component. We
   * can re-use it for the other vector components by rotating the
   * coordinates of the evaluation point.
   */
  const AnisotropicPolynomials<dim> polynomial_space;

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
PolynomialsABF<dim>::name() const
{
  return "ABF";
}


DEAL_II_NAMESPACE_CLOSE

#endif
