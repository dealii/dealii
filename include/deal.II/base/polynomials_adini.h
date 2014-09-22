// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2013 by the deal.II authors
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


#ifndef __deal2__polynomials_adini_h
#define __deal2__polynomials_adini_h

#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/table.h>

DEAL_II_NAMESPACE_OPEN

/**
 * The cubic polynomial space for the Adini element
 *
 * This space consists of the cubic space <i>P<sub>3</sub></i>
 * augmented by the functions <i>xy<sup>3</sup></i> and
 * <i>x<sup>3</sup>y</i>.
 *
 * The basis of the space is chosen to match the node functionals of
 * the Adini element.
 *
 * @todo This polynomial space is implemented in 2D only.
 *
 * @author Bärbel Janssen, 2007
 */

class PolynomialsAdini
{
public:
  /**
   * Constructor for
   * the polynomials of
   * the described space
   */
  PolynomialsAdini ();
  /**
   * Computes the value and the
   * first and second derivatives
   * of each polynomial at
   * <tt>unit_point</tt>.
   *
   * The size of the vectors must
   * either be equal 0 or equal
   * n(). In the first case,
   * the function will not compute
   * these values, i.e. you
   * indicate what you want to have
   * computed by resizing those
   * vectors which you want filled.
   *
   * If you need values or
   * derivatives of all polynomials
   * then use this function, rather
   * than using any of the
   * compute_value(),
   * compute_grad() or
   * compute_grad_grad()
   * functions, see below, in a
   * loop over all polynomials.
   */

  void compute (const Point<2> &unit_point,
                std::vector<double> &values,
                std::vector<Tensor<1,2> > &grads,
                std::vector< Tensor<2,2> > &grad_grads) const;

  /**
   * Computes the value of the
   * <tt>i</tt>th polynomial at
   * <tt>unit_point</tt>.
   *
   * Consider using compute() instead.
   */

  double compute_value (const unsigned int i,
                        const Point<2> &p) const;

  /**
   * Computes the gradient of the
   * <tt>i</tt>th polynomial at
   * <tt>unit_point</tt>.
   *
   * Consider using compute() instead.
   */

  Tensor<1,2> compute_grad (const unsigned int i,
                            const Point<2> &p) const;
  /**
   * Computes the second derivative
   * (grad_grad) of the <tt>i</tt>th
   * polynomial at
   * <tt>unit_point</tt>.
   *
   * Consider using compute() instead.
   */

  Tensor<2,2> compute_grad_grad (const unsigned int i, const Point<2> &p) const;
  Tensor<2,2> compute_grad_grad_2 (const unsigned int i, const Point<2> &p) const;

private:
  /**
   * Store the coefficients of the
   * polynominals in the order
   * $1,x,y,x^2,y^2,xy,x^3,y^3,xy^2,x^2y,x^3y,xy^3$
   */
  Table<2, double> coef;

  /**
   * Store the coefficients of the x-derivative
   * of the polynominals in the order
   * $1,x,y,x^2,y^2,xy,x^3,y^3,xy^2,x^2y,x^3y,xy^3$
   */

  Table<2, double> dx;
  /**
   * Store the coefficients of the y-derivative
   * of the polynominals in the order
   * $1,x,y,x^2,y^2,xy,x^3,y^3,xy^2,x^2y,x^3y,xy^3$
   */

  Table<2, double> dy;
  /**
   * Store the coefficients of the second x-derivative
   * of the polynominals in the order
   * $1,x,y,x^2,y^2,xy,x^3,y^3,xy^2,x^2y,x^3y,xy^3$
   */
  Table<2, double> dxx;
  /**
   * Store the coefficients of the second y-derivative
   * of the polynominals in the order
   * $1,x,y,x^2,y^2,xy,x^3,y^3,xy^2,x^2y,x^3y,xy^3$
   */
  Table<2, double> dyy;
  /**
   * Store the coefficients of the second mixed derivative
   * of the polynominals in the order
   * $1,x,y,x^2,y^2,xy,x^3,y^3,xy^2,x^2y,x^3y,xy^3$
   */
  Table<2, double> dxy;

};



DEAL_II_NAMESPACE_CLOSE

#endif
