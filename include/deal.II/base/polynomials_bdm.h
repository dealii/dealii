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

#ifndef __deal2__polynomials_BDM_h
#define __deal2__polynomials_BDM_h


#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/point.h>
#include <deal.II/base/polynomial.h>
#include <deal.II/base/polynomial_space.h>
#include <deal.II/base/table.h>
#include <deal.II/base/thread_management.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN

/**
 * This class implements the <i>H<sup>div</sup></i>-conforming,
 * vector-valued Brezzi-Douglas-Marini polynomials as described in the
 * book by Brezzi and Fortin.
 *
 * These polynomial spaces are based on the space
 * <i>P<sub>k</sub></i>, realized by a PolynomialSpace constructed
 * with Legendre polynomials. Since these shape functions are not
 * sufficient, additional functions are added. These are the following
 * vector valued polynomials:
 *
 * <dl>
 * <dt> In 2D:
 * <dd> The 2D-curl of the functions <i>x<sup>k+1</sup>y</i>
 * and <i>xy<sup>k+1</sup></i>.
 * <dt>In 3D:
 * <dd> For any <i>i=0,...,k</i> the curls of
 * <i>(0,0,xy<sup>i+1</sup>z<sup>k-i</sup>)</i>,
 * <i>(x<sup>k-i</sup>yz<sup>i+1</sup>,0,0)</i> and
 * <i>(0,x<sup>i+1</sup>y<sup>k-i</sup>z,0)</i>
 * </dl>
 *
 * @todo Second derivatives in 3D are missing.
 *
 * @ingroup Polynomials
 * @author Guido Kanschat
 * @date 2003, 2005, 2009
 */
template <int dim>
class PolynomialsBDM
{
public:
  /**
   * Constructor. Creates all basis
   * functions for BDM polynomials
   * of given degree.
   *
   * @arg k: the degree of the
   * BDM-space, which is the degree
   * of the largest complete
   * polynomial space
   * <i>P<sub>k</sub></i> contained
   * in the BDM-space.
   */
  PolynomialsBDM (const unsigned int k);

  /**
   * Computes the value and the
   * first and second derivatives
   * of each BDM
   * polynomial at @p unit_point.
   *
   * The size of the vectors must
   * either be zero or equal
   * <tt>n()</tt>.  In the
   * first case, the function will
   * not compute these values.
   *
   * If you need values or
   * derivatives of all tensor
   * product polynomials then use
   * this function, rather than
   * using any of the
   * <tt>compute_value</tt>,
   * <tt>compute_grad</tt> or
   * <tt>compute_grad_grad</tt>
   * functions, see below, in a
   * loop over all tensor product
   * polynomials.
   */
  void compute (const Point<dim>            &unit_point,
                std::vector<Tensor<1,dim> > &values,
                std::vector<Tensor<2,dim> > &grads,
                std::vector<Tensor<3,dim> > &grad_grads) const;

  /**
   * Returns the number of BDM polynomials.
   */
  unsigned int n () const;

  /**
   * Returns the degree of the BDM
   * space, which is one less than
   * the highest polynomial degree.
   */
  unsigned int degree () const;

  /**
   * Return the name of the space,
   * which is <tt>BDM</tt>.
   */
  std::string name () const;

  /**
   * Return the number of
   * polynomials in the space
   * <TT>BDM(degree)</tt> without
   * requiring to build an object
   * of PolynomialsBDM. This is
   * required by the FiniteElement
   * classes.
   */
  static unsigned int compute_n_pols(unsigned int degree);

private:
  /**
   * An object representing the
   * polynomial space used
   * here. The constructor fills
   * this with the monomial basis.
   */
  const PolynomialSpace<dim> polynomial_space;

  /**
   * Storage for monomials. In 2D,
   * this is just the polynomial of
   * order <i>k</i>. In 3D, we
   * need all polynomials from
   * degree zero to <i>k</i>.
   */
  std::vector<Polynomials::Polynomial<double> > monomials;

  /**
   * Number of BDM
   * polynomials.
   */
  unsigned int n_pols;

  /**
   * A mutex that guards the
   * following scratch arrays.
   */
  mutable Threads::Mutex mutex;

  /**
   * Auxiliary memory.
   */
  mutable std::vector<double> p_values;

  /**
   * Auxiliary memory.
   */
  mutable std::vector<Tensor<1,dim> > p_grads;

  /**
   * Auxiliary memory.
   */
  mutable std::vector<Tensor<2,dim> > p_grad_grads;
};


template <int dim>
inline unsigned int
PolynomialsBDM<dim>::n() const
{
  return n_pols;
}


template <int dim>
inline unsigned int
PolynomialsBDM<dim>::degree() const
{
  return polynomial_space.degree();
}


template <int dim>
inline std::string
PolynomialsBDM<dim>::name() const
{
  return "BDM";
}


DEAL_II_NAMESPACE_CLOSE

#endif
