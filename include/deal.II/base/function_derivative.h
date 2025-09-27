// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2000 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_function_derivative_h
#define dealii_function_derivative_h

#include <deal.II/base/config.h>

#include <deal.II/base/auto_derivative_function.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/function.h>

#include <vector>


DEAL_II_NAMESPACE_OPEN


/**
 * Derivative of a function object.  The value access functions of this class
 * return the directional derivative of a function with respect to a direction
 * provided on construction. If <tt>b</tt> is the vector, the derivative <tt>b
 * . grad f</tt> is computed. This derivative is evaluated directly, not by
 * computing the gradient of <tt>f</tt> and its scalar product with
 * <tt>b</tt>.
 *
 * The derivative is computed numerically, using one of the provided
 * difference formulas (see <tt>set_formula</tt> for available schemes).
 * Experimenting with <tt>h</tt> and the difference scheme may be necessary to
 * obtain sufficient results.
 *
 * @ingroup functions
 */
template <int dim>
class FunctionDerivative : public AutoDerivativeFunction<dim>
{
public:
  /**
   * Constructor. Provided are the functions to compute derivatives of, the
   * direction vector of the differentiation and the step size <tt>h</tt> of
   * the difference formula.
   */
  FunctionDerivative(const Function<dim> &f,
                     const Point<dim>    &direction,
                     const double         h = 1.e-6);

  /**
   * Constructor. Provided are the functions to compute derivatives of and the
   * direction vector of the differentiation in each quadrature point and the
   * difference step size.
   *
   * This is the constructor for a variable velocity field. Most probably, a
   * new object of <tt>FunctionDerivative</tt> has to be constructed for each
   * set of quadrature points.
   *
   * The number of quadrature point must still be the same, when values are
   * accessed.
   */
  FunctionDerivative(const Function<dim>           &f,
                     const std::vector<Point<dim>> &direction,
                     const double                   h = 1.e-6);

  /**
   * Choose the difference formula. This is set to the default in the
   * constructor.
   *
   * Formulas implemented right now are first order backward Euler
   * (<tt>UpwindEuler</tt>), second order symmetric Euler (<tt>Euler</tt>) and
   * a symmetric fourth order formula (<tt>FourthOrder</tt>).
   */
  void
  set_formula(typename AutoDerivativeFunction<dim>::DifferenceFormula formula =
                AutoDerivativeFunction<dim>::Euler);
  /**
   * Change the base step size of the difference formula
   */
  void
  set_h(const double h);

  virtual double
  value(const Point<dim> &p, const unsigned int component = 0) const override;

  virtual void
  vector_value(const Point<dim> &p, Vector<double> &value) const override;

  virtual void
  value_list(const std::vector<Point<dim>> &points,
             std::vector<double>           &values,
             const unsigned int             component = 0) const override;

  /**
   * Return an estimate for the memory consumption, in bytes, of this object.
   * This is not exact (but will usually be close) because calculating the
   * memory usage of trees (e.g., <tt>std::map</tt>) is difficult.
   */
  virtual std::size_t
  memory_consumption() const override;

private:
  /**
   * Function for differentiation.
   */
  const Function<dim> &f;

  /**
   * Step size of the difference formula.
   */
  double h;

  /**
   * Difference formula.
   */
  typename AutoDerivativeFunction<dim>::DifferenceFormula formula;

  /**
   * Helper object. Contains the increment vector for the formula.
   */
  std::vector<Tensor<1, dim>> incr;
};

DEAL_II_NAMESPACE_CLOSE

#endif
