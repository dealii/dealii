// ---------------------------------------------------------------------
//
// Copyright (C) 2019 - 2020 by the deal.II authors
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
 * the time step size during the simulation. This class provides the necessary
 * interface to be incorporated in any time-dependent simulation. As an
 * example, the usage of this class is demonstrated in step-21. This class
 * attempts to replace the usage of TimestepControl with a better and more
 * modern interface.
 *
 * This class provides a number of invariants that are guaranteed to be
 * true at all times.
 *
 * * The current simulation time is within the closed interval between the
 *   start time and the end time ($T_{\text{start}} \le t \le T_{\text{end}}$).
 * * Whenever time is incremented, the step size is positive ($dt > 0$).
 *   In other words, time advances in strictly ascending order
 *   ($m < n \Leftrightarrow t_m < t_n$).
 *
 * The model this class follows is that one sets a *desired* time step length
 * either through the constructor or using set_desired_next_step_size()
 * function. This step size will then be used in all following calls to the
 * advance_time() function, but may be adjusted slightly towards the end of the
 * simulation to ensure that the simulation time hits the end time exactly. The
 * adjustment is useful for the following reasons:
 *
 * Let's say that you loop over all of the time steps by using a `for` loop
 * @code
 *   for (DiscreteTime time(0., 1., 0.3);
 *        time.is_at_end() == false;
 *        time.advance_time())
 *   {
 *     // Insert simulation code here
 *   }
 * @endcode
 * or, if you like this better, the equivalent `while` loop:
 * @code
 *   DiscreteTime time(0., 1., 0.3);
 *   while (time.is_at_end() == false)
 *   {
 *     // Insert simulation code here
 *
 *     time.advance_time();
 *   }
 * @endcode
 *
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
 *        time.is_at_end() == false;
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
 *
 * ### Details of time-stepping
 *
 * Since time is marched forward in a discrete manner in our simulations, we
 * need to discuss how we increment time. During time stepping we enter two
 * separate alternating regimes in every step.
 *
 * * The **snapshot** stage (the **current** stage, the **consistent**
 *   stage): In this part of the algorithm, we are at $t = t_n$ and all
 *   quantities of the simulation (displacements, strains, temperatures, etc.)
 *   are up-to-date for $t = t_n$. In this stage, *current time* refers to
 *   $t_n$, *next time* refers to $t_{n+1}$, *previous time* refers to
 *   $t_{n-1}$. The other useful notation quantities are the *next* time step
 *   size $t_{n+1} - t_n$ and *previous* time step size $t_n - t_{n-1}$. In
 *   this stage, it is a perfect occasion to generate text output using print
 *   commands within the user's code. Additionally, post-processed outputs can
 *   be prepared here, which can then later be viewed by visualization programs
 *   such as `Tecplot`, `Paraview`, and `VisIt`. Additionally, during the
 *   snapshot stage, the code can assess the quality of the previous step and
 *   decide whether it wants to increase or decrease the time step size. The
 *   step size for the next time step can be modified here, by calling
 *   set_desired_next_step_size().
 * * The **update** stage (the **transition** stage, the **inconsistent**
 *   stage): In this section of the program, the internal state of the
 *   simulation is getting updated from $t_n$ to $t_{n+1}$. All of the
 *   variables need to be updated one by one, the step number is incremented,
 *   the time is incremented by $dt = t_{n+1} - t_n$, and time-integration
 *   algorithms are used to update the other simulation quantities. In the
 *   middle of this stage, some variables have been updated to $t_{n+1}$ but
 *   other variables still represent their value at $t_n$. Thus, we call this
 *   the inconsistent stage, requiring that no post-processing output related
 *   to the state variables take place within it. The state variables, namely
 *   those related to time, the solution field and any internal variables, are
 *   not synchronized and then get updated one by one. In general, the order of
 *   updating variables is arbitrary, but some care should be taken if there
 *   are interdependencies between them. For example, if some variable such as
 *   $x$ depends on the calculation of another variable such as $y$, then $y$
 *   must be updated before $x$ can be updated.
 *
 *   The question arises whether time should be incremented before updating
 *   state quantities. Multiple possibilities exist, depending on program and
 *   formulation requirements, and possibly the programmer's preferences:
 *   * Time is incremented *before* the rest of the updates. In this case, even
 *     though time is incremented to $t_{n+1}$, not all variables are updated
 *     yet. During this update phase, $dt$ equals the *previous* time step
 *     size. *Previous* means that it is referring to the $dt$ of the
 *     `advance_time()` command that was performed previously. In the
 *     following example code, we are assuming that `a` and `b` are two state
 *     variables that need to be updated in this time step.
 *     @code
 *       time.advance_time();
 *       new_a = update_a(a, b, time.get_previous_step_size());
 *       b = update_b(a, b, time.get_previous_step_size());
 *       a = new_a;
 *     @endcode
 *     Here, the code starts in a consistent state, but once advance_time()
 *     is called, the time variable, `a`, and `b` are no longer consistent
 *     with each other until after the last statement. At that point,
 *     the variables are all consistent again.
 *   * Time is incremented from $t_n$ to $t_{n+1}$ *after* all variables have
 *     already been updated for $t_{n+1}$. During the update stage, $dt$ is
 *     denoted as the *next* time step size. *Next* means that $dt$ of the
 *     step corresponds to the `advance_time()` command that will happen
 *     subsequently.
 *     @code
 *       new_a = update_a(a, b, time.get_next_step_size());
 *       b = update_b(a, b, time.get_next_step_size());
 *       a = new_a;
 *       time.advance_time();
 *     @endcode
 *   * Time is incremented in the middle of the other updates: In this case
 *     $dt$ would correspond to *next* or *previous* depending of whether it
 *     is used before or after the call to `advance_time()`.
 *     @code
 *       new_a = update_a(a, b, time.get_next_step_size());
 *       time.advance_time();
 *       b = update_b(a, b, time.get_previous_step_size());
 *       a = new_a;
 *     @endcode
 *
 * One thing to note is that, during the update phase, $dt$ is referred to
 * either **next** or **previous** time step size, depending on whether
 * advance_time() has been called yet. The notion of *current* time
 * step size is ill-defined. In fact, in the update stage the definition of
 * every variable depends on whether it has been updated yet or not, hence the
 * name **the inconsistent stage**.
 *
 * The following code snippet shows the code sections for the snapshot stage
 * and the update stage in the context of a complete time-dependent
 * simulation. This code follows the coding conventions incorporated in the
 * tutorial examples. Note that even though this example is written in the
 * format of a `for` loop, it can equivalently be written as a `while` or
 * `do while` loop (as shown in step-21).
 * @code
 * // pre-processing/setup stage {
 * make_grid();
 * setup_system();
 * for (DiscreteTime time(0., 1., 0.1);  // } end pre-processing/setup stage
 *      time.is_at_end() == false;
 *      time.advance_time())             // part of the update stage, runs at
 *                                       // the end of the loop body
 * {
 *   // snapshot stage {
 *   const double time_of_simulation = time.get_next_time();
 *   const double timestep_size      = time.get_next_step_size();
 *
 *   std::cout
 *     << "Timestep: " << time.get_step_number() << " -- "
 *     << "Solving for the solution at "
 *     << "t = " << time_of_simulation << " with "
 *     << "dt = " << timestep_size << "." << std::endl;
 *   // } end snapshot stage
 *
 *   // update stage {
 *   assemble_system(time_of_simulation, timestep_size);
 *   solve();
 *   update_solutions();
 *   // } end update stage
 *
 *   // snapshot stage {
 *   output_results(time_of_solution);
 *
 *   // propose a new timestep size if need be
 *   // time.set_desired_next_step_size(...);
 *   // } end snapshot stage
 * }
 * @endcode
 *
 * @author Reza Rastak, 2019
 */
class DiscreteTime
{
public:
  /**
   * Constructor.
   *
   * @param[in] start_time The time at the start of the simulation.
   *
   * @param[in] end_time The time at the end of the simulation.
   *
   * @param[in] desired_start_step_size A desired step size for incrementing
   * time for the first step. It is not guaranteed that this value will be
   * actually used as the size of the first step, as discussed in the
   * introduction.
   *
   * @pre @p desired_start_step_size must be non-negative.
   *
   * @note @p desired_start_step_size is an optional parameter. If it is not
   * provided or it is specified as zero, it indicates that the
   * desired size for the time step will be calculated at a different location
   * in the code. In this case, the created object cannot increment time until
   * the step size is changed by calling set_desired_next_step_size().
   */
  DiscreteTime(const double start_time,
               const double end_time,
               const double desired_start_step_size = 0.);

  /**
   * Return the current time.
   */
  double
  get_current_time() const;

  /**
   * Return the next time that we would reach if we were to advance the time
   * by one step.
   *
   * @note If the simulation is at the end time, this method returns the
   * end time.
   */
  double
  get_next_time() const;

  /**
   * Return the time we were at before `advance_time()` was called last time.
   *
   * @note If the simulation is at the start time, this method returns the
   * start time.
   */
  double
  get_previous_time() const;

  /**
   * Return the start time.
   */
  double
  get_start_time() const;

  /**
   * Return the end of the time interval.
   * The final time step ends exactly at this point. This exact floating-point
   * equality is very important because it allows us to equality-compare
   * current time with end time and decide whether we have reached the end of
   * the simulation.
   */
  double
  get_end_time() const;

  /**
   * Return whether no step has taken place yet.
   */
  bool
  is_at_start() const;

  /**
   * Return whether time has reached the end time.
   */
  bool
  is_at_end() const;

  /**
   * Return the size of the step from current time step to the
   * next. As discussed in the introduction to the class, this is the
   * *actual* time step, and may differ from the *desired* time step
   * set in the constructor or through the
   * set_desired_next_step_size() function.
   *
   * @note If the simulation is at the end time, this method returns zero.
   */
  double
  get_next_step_size() const;

  /**
   * Return the step size of the previous step.
   *
   * @note If the simulation is at the start time, this method returns zero.
   */
  double
  get_previous_step_size() const;

  /**
   * Return the number of times the simulation time has been incremented.
   * Return zero when the simulation is at the start time.
   */
  unsigned int
  get_step_number() const;

  /**
   * Set the *desired* value of the next time step size. By calling this
   * method, we are indicating the next time advance_time() is called, we
   * would like @p time_step_size to be used to advance the simulation time.
   * However, if the step is too large such that the next
   * simulation time exceeds the end time, the step size is truncated.
   * Additionally, if the step size is such that the next simulation time
   * approximates the end time (but falls just slightly short of it), the step
   * size is adjusted such that the next simulation time exactly matches the
   * end time.
   */
  void
  set_desired_next_step_size(const double time_step_size);

  /**
   * Advance the current time based on the value of the current step.
   * If you want to adjust the next time step size, call the method
   * set_desired_next_step_size() before calling this method.
   * If you call this function repeatedly, the time
   * is increased with the same step size until it reaches the end
   * time. See the documentation of set_desired_next_step_size() for
   * explanation of the rules for automatic adjustment of the step size.
   *
   * @pre Current time must be smaller than the end time. The object cannot
   * advance time if it is already at the end time. This rule is created to
   * avoid the creation of an infinite loop when advance_time() is called
   * inside a loop.
   *
   * @pre The time step size must be nonzero. If the step size is currently
   * zero, change it by calling set_desired_next_step_size() before calling
   * advance_time().
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
  double start_time;

  /**
   *The end of the time interval.
   */
  double end_time;

  /**
   * The current time.
   */
  double current_time;

  /**
   * The time at the next step.
   *
   * @note Internally, the next simulation time is stored instead of the
   * current step size. For example, when the method
   * set_desired_next_step_size() is called, it computes the appropriate next
   * simulation time and stores it. When advance_time() is called, the
   * current_time is replaced by next_time. This choice for the internal state
   * allows for simpler code and ensures than when we call advance_time() at
   * the last step, the floating-point value of the time exactly matches the
   * end time.
   */
  double next_time;

  /**
   * The previous time.
   */
  double previous_time;

  /**
   * The size of the first step.
   */
  double start_step_size;

  /**
   * The step number i.e. the number of times the simulation time ha been
   * incremented.
   */
  unsigned int step_number;
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



inline bool
DiscreteTime::is_at_start() const
{
  return step_number == 0;
}



inline bool
DiscreteTime::is_at_end() const
{
  return current_time == end_time;
}



inline double
DiscreteTime::get_next_step_size() const
{
  return next_time - current_time;
}



inline double
DiscreteTime::get_previous_step_size() const
{
  return current_time - previous_time;
}



inline double
DiscreteTime::get_current_time() const
{
  return current_time;
}



inline double
DiscreteTime::get_next_time() const
{
  return next_time;
}



inline double
DiscreteTime::get_previous_time() const
{
  return previous_time;
}



inline unsigned int
DiscreteTime::get_step_number() const
{
  return step_number;
}


DEAL_II_NAMESPACE_CLOSE

#endif
