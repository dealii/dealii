// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 1999 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_function_time_h
#define dealii_function_time_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>

DEAL_II_NAMESPACE_OPEN
/**
 * Support for time dependent functions. The library was also designed for
 * time dependent problems. For this purpose, the function objects also
 * contain a field which stores the time, as well as functions manipulating
 * them. Time independent problems should not access or even abuse them for
 * other purposes, but since one normally does not create thousands of
 * function objects, the gain in generality weighs out the fact that we need
 * not store the time value for not time dependent problems. The second
 * advantage is that the derived standard classes like <tt>ZeroFunction</tt>,
 * <tt>ConstantFunction</tt> etc also work for time dependent problems.
 *
 * Access to the time goes through the following functions:
 * <ul>
 *  <li> <tt>get_time</tt>: return the present value of the time variable.
 *  <li> <tt>set_time</tt>: set the time value to a specific value.
 *  <li> <tt>advance_time</tt>: increase the time by a certain time step.
 * </ul>
 * The latter two functions are virtual, so that derived classes can perform
 * computations which need only be done once for every new time. For example,
 * if a time dependent function had a factor <tt>sin(t)</tt>, then it may be a
 * reasonable choice to calculate this factor in a derived version of
 * set_time(), store it in a member variable and use that one rather than
 * computing it every time <tt>value()</tt>, <tt>value_list</tt> or one of the
 * other functions of class Function is called.
 *
 * By default, the advance_time() function calls the set_time() function with
 * the new time, so it is sufficient in most cases to overload only set_time()
 * for computations as sketched out above.
 *
 * The constructor of this class takes an initial value for the time variable,
 * which defaults to zero. Because a default value is given, none of the
 * derived classes needs to take an initial value for the time variable if not
 * needed.
 *
 * @tparam Number The data type in which time values are to be stored. This
 * will, in almost all cases, simply be the default @p double, but there are
 * cases where one may want to store the time in a different (and always
 * scalar) type. An example would be an interval type that can store a value
 * as well as its uncertainty. Another example would be a type that allows for
 * Automatic Differentiation (see, for example, the Sacado type used in
 * step-33) and thereby can generate analytic (temporal) derivatives of a
 * function.
 *
 *
 * @ingroup functions
 */
template <typename Number = double>
class FunctionTime
{
public:
  /**
   * Constructor. May take an initial value for the time variable, which
   * defaults to zero.
   */
  FunctionTime(const Number initial_time = Number(0.0));

  /**
   * Virtual destructor.
   */
  virtual ~FunctionTime() = default;

  /**
   * Return the value of the time variable.
   */
  Number
  get_time() const;

  /**
   * Set the time to <tt>new_time</tt>, overwriting the old value.
   */
  virtual void
  set_time(const Number new_time);

  /**
   * Advance the time by the given time step <tt>delta_t</tt>.
   */
  virtual void
  advance_time(const Number delta_t);

  /**
   * The type this class is initialized with and that is used to represent time.
   */
  using time_type = Number;

private:
  /**
   * Store the present time.
   */
  Number time;
};



/*----------------------------- Inline functions ----------------------------*/

#ifndef DOXYGEN

template <typename Number>
inline Number
FunctionTime<Number>::get_time() const
{
  return time;
}

#endif
DEAL_II_NAMESPACE_CLOSE

#endif
