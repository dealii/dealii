// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2004 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_polynomials_raviart_thomas_h
#define dealii_polynomials_raviart_thomas_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/point.h>
#include <deal.II/base/polynomial.h>
#include <deal.II/base/polynomial_space.h>
#include <deal.II/base/polynomials_vector_anisotropic.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/tensor_product_polynomials.h>

#include <mutex>
#include <vector>


DEAL_II_NAMESPACE_OPEN

/**
 * This class implements the <i>H<sup>div</sup></i>-conforming,
 * Raviart-Thomas polynomials as described in the book by Brezzi and Fortin.
 * Most of the functionality comes from the vector-valued anisotropic
 * polynomials class PolynomialsVectorAnisotropic.
 *
 * The Raviart-Thomas polynomials are constructed such that the divergence is
 * in the tensor product polynomial space <i>Q<sub>k</sub></i>. Therefore, the
 * polynomial order of each component must be one order higher in the
 * corresponding direction, yielding the polynomial spaces
 * <i>(Q<sub>k+1,k</sub>, Q<sub>k,k+1</sub>)</i> and <i>(Q<sub>k+1,k,k</sub>,
 * Q<sub>k,k+1,k</sub>, Q<sub>k,k,k+1</sub>)</i> in 2d and 3d, resp.
 *
 * @ingroup Polynomials
 */
template <int dim>
class PolynomialsRaviartThomas : public PolynomialsVectorAnisotropic<dim>
{
public:
  /**
   * Constructor, using the common Raviart-Thomas space of degree `k + 1` in
   * normal direction and `k` in the tangential directions.
   *
   * @arg k: the degree of the Raviart-Thomas-space, which is the degree of
   * the largest tensor product polynomial space <i>Q<sub>k</sub></i>
   * contains.
   */
  PolynomialsRaviartThomas(const unsigned int k);

  /**
   * Variant of the n_polynomials() function taking only a single argument
   * `degree`, assuming `degree + 1` in the normal direction and `degree` in
   * the tangential directions.
   */
  static unsigned int
  n_polynomials(const unsigned int degree);

  // Make respective two-argument method from base class available
  using PolynomialsVectorAnisotropic<dim>::n_polynomials;

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
};


DEAL_II_NAMESPACE_CLOSE

#endif
