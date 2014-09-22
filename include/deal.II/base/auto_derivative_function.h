// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2013 by the deal.II authors
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

#ifndef __deal2__auto_derivative_function_h
#define __deal2__auto_derivative_function_h


#include <deal.II/base/exceptions.h>
#include <deal.II/base/function.h>

DEAL_II_NAMESPACE_OPEN

/**
 * This class automatically computes the gradient of a function by
 * employing numerical difference quotients. This only, if the user
 * function does not provide the gradient function himself.
 *
 * The following example of an user defined function overloads and
 * implements only the value() function but not the gradient()
 * function. If the gradient() function is invoked then the gradient
 * function implemented by the AutoDerivativeFunction is called,
 * where the latter function imployes numerical difference quotients.
 *
 * @code
 * class UserFunction: public AutoDerivativeFunction
 * {               // access to one component at one point
 *   double value (const Point<dim> &p, const
 *                 unsigned int component = 0) const
 *          { // Implementation ....  };
 * } user_function;
 *
 *            // gradient by employing difference quotients.
 * Tensor<1,dim> grad=user_function.gradient(some_point);
 * @endcode
 *
 * If the user overloads and implements also the gradient function,
 * then, of course, the users gradient function is called.
 *
 * Note, that the usage of the value() and gradient() functions
 * explained above, also applies to the value_list() and
 * gradient_list() functions as well as to the vector valued
 * versions of these functions, see e.g. vector_value(),
 * vector_gradient(), vector_value_list() and
 * vector_gradient_list().
 *
 * The gradient() and gradient_list() functions make use of the
 * Function::value() function. The vector_gradient() and
 * vector_gradient_list() make use of the Function::vector_value()
 * function. Make sure that the user defined function implements the
 * value() function and the vector_value() function, respectively.
 *
 * Furthermore note, that an object of this class does <b>not</b> represent
 * the derivative of a function, like FunctionDerivative, that
 * gives a directional derivate by calling the value() function. In
 * fact, this class (the AutoDerivativeFunction class) can
 * substitute the Function class as base class for user defined
 * classes. This class implements the gradient() functions for
 * automatic computation of numerical difference quotients and serves
 * as intermediate class between the base Function class and the
 * user defined function class.
 *
 * @ingroup functions
 * @author Ralf Hartmann, 2001
 */
template <int dim>
class AutoDerivativeFunction : public Function<dim>
{
public:

  /**
   * Names of difference formulas.
   */
  enum DifferenceFormula
  {
    /**
     * The symmetric Euler
     * formula of second order:
     * @f[
     * u'(t) \approx
     * \frac{u(t+h) -
     * u(t-h)}{2h}.
     * @f]
     */
    Euler,
    /**
     * The upwind Euler
     * formula of first order:
     * @f[
     * u'(t) \approx
     * \frac{u(t) -
     * u(t-h)}{h}.
     * @f]
     */
    UpwindEuler,
    /**
     * The fourth order scheme
     * @f[
     * u'(t) \approx
     * \frac{u(t-2h) - 8u(t-h)
     * +  8u(t+h) - u(t+2h)}{12h}.
     * @f]
     */
    FourthOrder
  };

  /**
   * Constructor. Takes the
   * difference step size
   * <tt>h</tt>. It's within the user's
   * responsibility to choose an
   * appropriate value here. <tt>h</tt>
   * should be chosen taking into
   * account the absolute value as
   * well as the amount of local
   * variation of the function.
   * Setting <tt>h=1e-6</tt> might be a
   * good choice for functions with
   * an absolute value of about 1,
   * that furthermore does not vary
   * to much.
   *
   * <tt>h</tt> can be changed later
   * using the set_h() function.
   *
   * Sets DifferenceFormula
   * <tt>formula</tt> to the default
   * <tt>Euler</tt> formula of the
   * set_formula()
   * function. Change this preset
   * formula by calling the
   * set_formula() function.
   */
  AutoDerivativeFunction (const double h,
                          const unsigned int n_components = 1,
                          const double       initial_time = 0.0);

  /**
   * Virtual destructor; absolutely
   * necessary in this case.
   */
  virtual ~AutoDerivativeFunction ();

  /**
   * Choose the difference formula.
   * See the enum #DifferenceFormula
   * for available choices.
   */
  void set_formula (const DifferenceFormula formula = Euler);

  /**
   * Takes the difference step size
   * <tt>h</tt>. It's within the user's
   * responsibility to choose an
   * appropriate value here. <tt>h</tt>
   * should be chosen taking into
   * account the absolute value of
   * as well as the amount of local
   * variation of the function.
   * Setting <tt>h=1e-6</tt> might be a
   * good choice for functions with
   * an absolute value of about 1,
   * that furthermore does not vary
   * to much.
   */
  void set_h (const double h);

  /**
   * Return the gradient of the
   * specified component of the
   * function at the given point.
   *
   * Computes numerical difference
   * quotients using the preset
   * #DifferenceFormula.
   */
  virtual Tensor<1,dim> gradient (const Point<dim>   &p,
                                  const unsigned int  component = 0) const;

  /**
   * Return the gradient of all
   * components of the
   * function at the given point.
   *
   * Computes numerical difference
   * quotients using the preset
   * #DifferenceFormula.
   */
  virtual void vector_gradient (const Point<dim>            &p,
                                std::vector<Tensor<1,dim> > &gradients) const;

  /**
   * Set <tt>gradients</tt> to the
   * gradients of the specified
   * component of the function at
   * the <tt>points</tt>.  It is assumed
   * that <tt>gradients</tt> already has the
   * right size, i.e.  the same
   * size as the <tt>points</tt> array.
   *
   * Computes numerical difference
   * quotients using the preset
   * #DifferenceFormula.
   */
  virtual void gradient_list (const std::vector<Point<dim> > &points,
                              std::vector<Tensor<1,dim> >    &gradients,
                              const unsigned int              component = 0) const;

  /**
   * Set <tt>gradients</tt> to the gradients of
   * the function at the <tt>points</tt>,
   * for all components.
   * It is assumed that <tt>gradients</tt>
   * already has the right size, i.e.
   * the same size as the <tt>points</tt> array.
   *
   * The outer loop over
   * <tt>gradients</tt> is over the points
   * in the list, the inner loop
   * over the different components
   * of the function.
   *
   * Computes numerical difference
   * quotients using the preset
   * #DifferenceFormula.
   */
  virtual void vector_gradient_list (const std::vector<Point<dim> > &points,
                                     std::vector<std::vector<Tensor<1,dim> > > &gradients) const;

  /**
   * Returns a
   * #DifferenceFormula of the
   * order <tt>ord</tt> at minimum.
   */
  static
  DifferenceFormula
  get_formula_of_order (const unsigned int ord);

  /**
   * Exception.
   */
  DeclException0(ExcInvalidFormula);

private:

  /**
   * Step size of the difference
   * formula. Set by the set_h()
   * function.
   */
  double h;

  /**
   * Includes the unit vectors
   * scaled by <tt>h</tt>.
   */
  std::vector<Tensor<1,dim> > ht;

  /**
   * Difference formula. Set by the
   * set_formula() function.
   */
  DifferenceFormula formula;
};


DEAL_II_NAMESPACE_CLOSE

#endif
