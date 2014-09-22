// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2013 by the deal.II authors
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

#ifndef __deal2__function_time_h
#define __deal2__function_time_h


#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>

DEAL_II_NAMESPACE_OPEN
/**
 *  Support for time dependent functions.
 *  The library was also designed for time dependent problems. For this
 *  purpose, the function objects also contain a field which stores the
 *  time, as well as functions manipulating them. Time independent problems
 *  should not access or even abuse them for other purposes, but since one
 *  normally does not create thousands of function objects, the gain in
 *  generality weighs out the fact that we need not store the time value
 *  for not time dependent problems. The second advantage is that the derived
 *  standard classes like <tt>ZeroFunction</tt>, <tt>ConstantFunction</tt> etc also work
 *  for time dependent problems.
 *
 *  Access to the time goes through the following functions:
 *  @verbatim
 *  <li> <tt>get_time</tt>: return the present value of the time variable.
 *  <li> <tt>set_time</tt>: set the time value to a specific value.
 *  <li> <tt>advance_time</tt>: increase the time by a certain time step.
 *  @endverbatim
 *  The latter two functions are virtual, so that derived classes can
 *  perform computations which need only be done once for every new time.
 *  For example, if a time dependent function had a factor <tt>sin(t)</tt>, then
 *  it may be a reasonable choice to calculate this factor in a derived
 *  version of <tt>set_time</tt>, store it in a member variable and use that one
 *  rather than computing it every time <tt>operator()</tt>, <tt>value_list</tt> or one
 *  of the other functions is called.
 *
 *  By default, the <tt>advance_time</tt> function calls the <tt>set_time</tt> function
 *  with the new time, so it is sufficient in most cases to overload only
 *  <tt>set_time</tt> for computations as sketched out above.
 *
 *  The constructor of this class takes an initial value for the time
 *  variable, which defaults to zero. Because a default value is given,
 *  none of the derived classes needs to take an initial value for the
 *  time variable if not needed.
 *
 *  Once again the warning: do not use the <tt>time</tt> variable for any other
 *  purpose than the intended one! This will inevitably lead to confusion.
 *
 *
 *  @ingroup functions
 *  @author Wolfgang Bangerth, Guido Kanschat, 1998, 1999
 */
template <typename Number=double>
class FunctionTime
{
public:
  /**
   * Constructor. May take an initial vakue
   * for the time variable, which defaults
   * to zero.
   */
  FunctionTime (const Number initial_time = Number(0.0));

  /**
   * Virtual destructor.
   */
  virtual ~FunctionTime();

  /**
   * Return the value of the time variable/
   */
  Number get_time () const;

  /**
   * Set the time to <tt>new_time</tt>, overwriting
   * the old value.
   */
  virtual void set_time (const Number new_time);

  /**
   * Advance the time by the given
   * time step <tt>delta_t</tt>.
   */
  virtual void advance_time (const Number delta_t);

private:
  /**
   * Store the present time.
   */
  Number time;
};



/*------------------------------ Inline functions ------------------------------*/

#ifndef DOXYGEN

template<typename Number>
inline Number
FunctionTime<Number>::get_time () const
{
  return time;
}

#endif
DEAL_II_NAMESPACE_CLOSE

#endif
