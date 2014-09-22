// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#ifndef __deal2__polynomials_P_h
#define __deal2__polynomials_P_h


#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/point.h>
#include <deal.II/base/polynomial.h>
#include <deal.II/base/polynomial_space.h>
#include <deal.II/base/table.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN
/**
 * @addtogroup Polynomials
 * @{
 */

/**
 * This class implements the polynomial space of degree <tt>p</tt>
 * based on the monomials ${1,x,x^2,...}$. I.e. in <tt>d</tt>
 * dimensions it constructs all polynomials of the form $\prod_{i=1}^d
 * x_i^{n_i}$, where $\sum_i n_i\leq p$. The base polynomials are
 * given a specific ordering, e.g. in 2 dimensions:
 * ${1,x,y,xy,x^2,y^2,x^2y,xy^2,x^3,y^3,...}$. The ordering of the
 * monomials in $P_k1$ matches the ordering of the monomials in $P_k2$
 * for $k2>k1$.
 *
 * @author Ralf Hartmann, 2004
 */
template <int dim>
class PolynomialsP: public PolynomialSpace<dim>
{
public:
  /**
   * Access to the dimension of
   * this object, for checking and
   * automatic setting of dimension
   * in other classes.
   */
  static const unsigned int dimension = dim;

  /**
   * Constructor. Creates all basis
   * functions of $P_p$.
   * @arg p: the degree of the
   * polynomial space
   */
  PolynomialsP (const unsigned int p);

  /**
   * Returns the degree <tt>p</tt>
   * of the polynomial space
   * <tt>P_p</tt>.
   *
   * Note, that this number is
   * <tt>PolynomialSpace::degree()-1</tt>,
   * compare definition in
   * PolynomialSpace.
   */
  unsigned int degree() const;

  /**
   * For the <tt>n</tt>th
   * polynomial $p_n(x,y,z)=x^i y^j
   * z^k$ this function gives the
   * degrees i,j,k in the x,y,z
   * directions.
   */
  void directional_degrees(unsigned int n,
                           unsigned int (&degrees)[dim]) const;

private:

  /**
   * Fills the <tt>index_map</tt>.
   */
  void create_polynomial_ordering(std::vector<unsigned int> &index_map) const;

  /**
   * Degree <tt>p</tt> of the
   * polynomial space $P_p$,
   * i.e. the number <tt>p</tt>
   * which was given to the
   * constructor.
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
inline void
PolynomialsP<dim>::directional_degrees(unsigned int n,
                                       unsigned int (&degrees)[dim]) const
{
  this->compute_index(n,degrees);
}

DEAL_II_NAMESPACE_CLOSE

#endif
