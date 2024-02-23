// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2004 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_polynomials_P_h
#define dealii_polynomials_P_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/point.h>
#include <deal.II/base/polynomial.h>
#include <deal.II/base/polynomial_space.h>
#include <deal.II/base/tensor.h>

#include <memory>
#include <vector>

DEAL_II_NAMESPACE_OPEN
/**
 * @addtogroup Polynomials
 * @{
 */

/**
 * This class implements the polynomial space of degree <tt>p</tt> based on
 * the monomials ${1,x,x^2,...}$. I.e. in <tt>d</tt> dimensions it constructs
 * all polynomials of the form $\prod_{i=1}^d x_i^{n_i}$, where $\sum_i
 * n_i\leq p$. The base polynomials are given a specific ordering, e.g. in 2
 * dimensions: ${1,x,y,xy,x^2,y^2,x^2y,xy^2,x^3,y^3,...}$. The ordering of the
 * monomials in $P_k1$ matches the ordering of the monomials in $P_k2$ for
 * $k2>k1$.
 */
template <int dim>
class PolynomialsP : public PolynomialSpace<dim>
{
public:
  /**
   * Access to the dimension of this object, for checking and automatic
   * setting of dimension in other classes.
   */
  static constexpr unsigned int dimension = dim;

  /**
   * Constructor. Creates all basis functions of $P_p$. @arg p: the degree of
   * the polynomial space
   */
  PolynomialsP(const unsigned int p);

  /**
   * Return the degree <tt>p</tt> of the polynomial space <tt>P_p</tt>.
   *
   * Note, that this number is <tt>PolynomialSpace::degree()-1</tt>, compare
   * definition in PolynomialSpace.
   */
  virtual unsigned int
  degree() const override;

  /**
   * For the <tt>n</tt>th polynomial $p_n(x,y,z)=x^i y^j z^k$ this function
   * gives the degrees i,j,k in the x,y,z directions.
   *
   * In 1d and 2d, obviously only i and i,j are returned.
   */
  std::array<unsigned int, dim>
  directional_degrees(unsigned int n) const;

  std::unique_ptr<ScalarPolynomialsBase<dim>>
  clone() const override
  {
    return std::make_unique<PolynomialsP<dim>>(*this);
  }

private:
  /**
   * Fills the <tt>index_map</tt>.
   */
  void
  create_polynomial_ordering(std::vector<unsigned int> &index_map) const;

  /**
   * Degree <tt>p</tt> of the polynomial space $P_p$, i.e. the number
   * <tt>p</tt> which was given to the constructor.
   */
  const unsigned int p;
};

/** @} */

template <int dim>
inline unsigned int
PolynomialsP<dim>::degree() const
{
  return p;
}


template <int dim>
inline std::array<unsigned int, dim>
PolynomialsP<dim>::directional_degrees(unsigned int n) const
{
  return this->compute_index(n);
}

DEAL_II_NAMESPACE_CLOSE

#endif
