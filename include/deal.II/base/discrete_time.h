// ---------------------------------------------------------------------
//
// Copyright (C) 2019 by the deal.II authors
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

#ifndef dealii_discrete_time_h
#define dealii_discrete_time_h

#include <deal.II/base/config.h>

DEAL_II_NAMESPACE_OPEN

/**
 * This class provides a means to keep track of the simulation time in a
 * time-dependent simulation. It manages stepping forward from a start time
 * $T_{\text{start}}$ to an end time $T_{\text{end}}$. It also allows adjusting
 * the time step size during the simulation. It is guaranteed that at all times
 * the current simulation time is in the closed interval between the start time
 * and the end time.
 *
 * The model this class follows is that one sets a *desired* time step length
 * either through the constructor or using set_desired_next_step_size()
 * function. This step size will then be used in all following calls to the
 * advance_time() function, but may be adjusted slightly towards the end of the
 * simulation to ensure that the simulation time hits the end time exactly. The
 * adjustment is useful for the following reasons:
 *
 * Let's say that you loop over all of the time steps by using a for loop
 * @code
 *   for (DiscreteTime time(0., 1., 0.3);
 *        time.get_current_time() != time.get_end_time();
 *        time.advance_time())
 *   {
 *     // Insert simulation code here
 *   }
 * @endcode
 * In the above example the time starts at $T_{\text{start}} = 0$ until
 * $T_{\text{end}}=1$. Assuming the time step $dt = 0.3$ is not modified inside
 * the loop, the time is advanced from $t = 0$ to $t = 0.3$, $t = 0.6$, $t =
 * 0.9$ and finally it reaches the end time at $t = 1.0$. Here, the final step
 * size needs to be reduced from its desired value of 0.3 to $dt = 0.1$ in order
 * to ensure that we finish the simulation exactly at the specified end time. In
 * fact, you should assume that not only the last time step length may be
 * adjusted, but also previously ones -- for example, this class may take the
 * liberty to spread the decrease in time step size out over several time steps
 * and increment time from $t=0$, to $0.3$, $0.6$, $0.8$, and finally
 * $t=T_{\text{end}}=1$ to avoid too large a change in time step size from one
 * step to another.
 *
 * The other situation in which the time step needs to be adjusted (this time to
 * slightly larger values) is if a time increment falls just short of the final
 * time. Imagine, for example, a similar situation as above, but with different
 * end time:
 * @code
 *   for (DiscreteTime time(0., 1.21, 0.3);
 *        time.get_current_time() != time.get_end_time();
 *        time.advance_time())
 *   {
 *     // Insert simulation code here
 *   }
 * @endcode
 * Here, the time step from $t=0.9$ to $t=1.2$ falls just short of the final
 * time $T_{\text{end}}=1.21$. Instead of following up with a very small step of
 * length $dt=0.01$, the class stretches the last time step (or last time steps)
 * slightly to reach the desired end time.
 *
 * The examples above make clear that the time step size given to this class is
 * only a *desired* step size. You can query the actual time step size using the
 * get_next_step_size() function.
 *
 * @author Reza Rastak, 2019
 */
class DiscreteTime
{
public:
  /**
   * Constructor
   */
  DiscreteTime(const double start_time,
               const double end_time,
               const double start_step_size);

  /**
   * Return the current time.
   */
  double
  get_current_time() const;

  /**
   * Return the start time.
   */
  double
  get_start_time() const;

  /**
   * Return the end of the time interval.
   * The final time step ends exactly at this point. This exact floating-point
   * equality is very important because it allows us to use the expression
   * <code>time.get_current_time() != time.get_end_time()</code> as the
   * conditional statement in a for loop to check if the end time is reached.
   */
  double
  get_end_time() const;

  /**
   * Return the size of the step from current time step to the
   * next. As discussed in the introduction to the class, this is the
   * *actual* time step, and may differ from the *desired* time step
   * set in the constructor or through the
   * set_desired_next_step_size() function.
   */
  double
  get_next_step_size() const;

  /**
   * Set the value of the next time step size. The next time advance_time()
   * is called, the newly set @p time_step_size will be used to advance
   * the simulation time. However, if the step is too large such that the next
   * simulation time exceeds the end time, the step size is truncated.
   * Additionally, if the step size is such that the next simulation time
   * approximates the end time (but falls just slightly short of it), the step
   * size is adjusted such that the next simulation time exactly matches the
   * end time.
   */
  void
  set_next_step_size(const double time_step_size);

  /**
   * Advance the current time based on the value of the current step.
   * If you want to adjust the next time step size, call the method
   * set_next_step_size() before calling this method.
   * If you call this function repeatedly, the time
   * is increased with the same step size until it reaches the end
   * time. See the documentation of set_next_step_size() for explanation
   * of the rules for automatic adjustment of the step size.
   *
   * @pre Current time must be smaller than the end time. The object cannot
   * advance time if it is already at the end time. This rule is created to
   * avoid the creation of an infinite loop when advance_time() is called
   * inside a loop.
   */
  void
  advance_time();

  /**
   * Set the current time equal to start time and set the step size to the
   * initial step size.
   */
  void
  restart();

private:
  /**
   * The beginning of the time interval.
   */
  const double start_time;

  /**
   *The end of the time interval.
   */
  const double end_time;

  /**
   * The size of the first step.
   */
  const double start_step_size;

  /**
   * The current time.
   */
  double current_time;

  /**
   * The time at the next step.
   *
   * @note Internally, the next simulation time is stored instead of the
   * current step size. For example, when the method set_next_step_size()
   * is called, it computes the appropriate next simulation time and stores
   * it. When advance_time() is called, the current_time is replaced by
   * next_time. This choice for the internal state allows for simpler code
   * and ensures than when we call advance_time() at the last step, the
   * floating-point value of the time exactly matches the end time.
   */
  double next_time;
};


/*---------------------- Inline functions ------------------------------*/


inline double
DiscreteTime::get_start_time() const
{
  return start_time;
}



inline double
DiscreteTime::get_end_time() const
{
  return end_time;
}



inline double
DiscreteTime::get_next_step_size() const
{
  return next_time - current_time;
}



inline double
DiscreteTime::get_current_time() const
{
  return current_time;
}


DEAL_II_NAMESPACE_CLOSE

#endif
