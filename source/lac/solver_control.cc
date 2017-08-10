// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2016 by the deal.II authors
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

#include <deal.II/base/logstream.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/signaling_nan.h>
#include <deal.II/lac/solver_control.h>

#include <cmath>
#include <sstream>

DEAL_II_NAMESPACE_OPEN

/*----------------------- SolverControl ---------------------------------*/


SolverControl::SolverControl (const unsigned int maxiter,
                              const double tolerance,
                              const bool m_log_history,
                              const bool m_log_result)
  :
  maxsteps(maxiter),
  tol(tolerance),
  lcheck(failure),
  initial_val(numbers::signaling_nan<double>()),
  lvalue(numbers::signaling_nan<double>()),
  lstep(0),
  check_failure(false),
  relative_failure_residual(0),
  failure_residual(0),
  history_data_enabled(false)
{}



SolverControl::~SolverControl()
{}



SolverControl::State
SolverControl::check (const unsigned int step,
                      const double check_value)
{
  // if this is the first time we
  // come here, then store the
  // residual for later comparisons
  if (step==0)
    {
      initial_val = check_value;
    }

  lstep  = step;
  lvalue = check_value;

  if (step==0)
    {
      if (check_failure)
        failure_residual=relative_failure_residual*check_value;
    }

  if (history_data_enabled)
    history_data.push_back(check_value);

  if (check_value <= tol)
    {
      lcheck = success;
      return success;
    }

  if ((step >= maxsteps) ||
      numbers::is_nan(check_value) ||
      (check_failure && (check_value > failure_residual))
     )
    {
      lcheck = failure;
      return failure;
    }

  lcheck = iterate;
  return iterate;
}



SolverControl::State
SolverControl::last_check() const
{
  return lcheck;
}


double
SolverControl::initial_value() const
{
  return initial_val;
}


double
SolverControl::last_value() const
{
  return lvalue;
}


unsigned int
SolverControl::last_step() const
{
  return lstep;
}


void
SolverControl::enable_history_data ()
{
  history_data_enabled = true;
}



const std::vector<double> &SolverControl::get_history_data() const
{
  Assert (history_data_enabled, ExcHistoryDataRequired());
  Assert (history_data.size() > 0,
          ExcMessage("The SolverControl object was asked for the solver history "
                     "data, but there is no data. Possibly you requested the data before the "
                     "solver was run."));

  return history_data;
}



double
SolverControl::average_reduction() const
{
  if (lstep == 0)
    return 0.;

  Assert (history_data_enabled, ExcHistoryDataRequired());
  Assert (history_data.size() > lstep, ExcInternalError());
  Assert (history_data[0] > 0., ExcInternalError());
  Assert (history_data[lstep] > 0., ExcInternalError());

  return std::pow(history_data[lstep]/history_data[0], 1./lstep);
}



double
SolverControl::step_reduction(unsigned int step) const
{
  Assert (history_data_enabled, ExcHistoryDataRequired());
  Assert (history_data.size() > lstep, ExcInternalError());
  Assert (step <=lstep, ExcIndexRange(step,1,lstep+1));
  Assert (step>0, ExcIndexRange(step,1,lstep+1));

  return history_data[step]/history_data[step-1];
}


double
SolverControl::final_reduction() const
{
  return step_reduction(lstep);
}


void
SolverControl::declare_parameters (ParameterHandler &param)
{
  param.declare_entry ("Max steps", "100", Patterns::Integer());
  param.declare_entry ("Tolerance", "1.e-10", Patterns::Double());
}


void SolverControl::parse_parameters (ParameterHandler &param)
{
  set_max_steps (param.get_integer("Max steps"));
  set_tolerance (param.get_double("Tolerance"));
}

/*----------------------- ReductionControl ---------------------------------*/


ReductionControl::ReductionControl(const unsigned int n,
                                   const double tol,
                                   const double red,
                                   const bool m_log_history,
                                   const bool m_log_result)
  :
  SolverControl (n, tol, m_log_history, m_log_result),
  reduce(red),
  reduced_tol(numbers::signaling_nan<double>())
{}


ReductionControl::ReductionControl (const SolverControl &c)
  :
  SolverControl(c),
  reduce(numbers::signaling_nan<double>()),
  reduced_tol(numbers::signaling_nan<double>())
{
  set_reduction(0.);
}


ReductionControl &
ReductionControl::operator= (const SolverControl &c)
{
  SolverControl::operator=(c);
  set_reduction(0.);
  return *this;
}


ReductionControl::~ReductionControl()
{}


SolverControl::State
ReductionControl::check (const unsigned int step,
                         const double check_value)
{
  // if this is the first time we
  // come here, then store the
  // residual for later comparisons
  if (step==0)
    {
      initial_val = check_value;
      reduced_tol = check_value * reduce;
    }

  // check whether desired reduction
  // has been achieved. also check
  // for equality in case initial
  // residual already was zero
  if (check_value <= reduced_tol)
    {
      lstep  = step;
      lvalue = check_value;

      lcheck = success;
      return success;
    }
  else
    return SolverControl::check(step, check_value);
}



void
ReductionControl::declare_parameters (ParameterHandler &param)
{
  SolverControl::declare_parameters (param);
  param.declare_entry("Reduction", "1.e-2", Patterns::Double());
}


void
ReductionControl::parse_parameters (ParameterHandler &param)
{
  SolverControl::parse_parameters (param);
  set_reduction (param.get_double("Reduction"));
}


/*---------------------- IterationNumberControl -----------------------------*/


IterationNumberControl::IterationNumberControl(const unsigned int n,
                                               const double       tolerance,
                                               const bool m_log_history,
                                               const bool m_log_result)
  :
  SolverControl (n, tolerance, m_log_history, m_log_result) {}


IterationNumberControl::~IterationNumberControl()
{}


SolverControl::State
IterationNumberControl::check (const unsigned int step,
                               const double check_value)
{
  // check whether the given number of iterations was reached, and return
  // success in that case. Otherwise, go on to the check of the base class.
  if (step >= this->maxsteps)
    {
      lstep  = step;
      lvalue = check_value;

      lcheck = success;
      return success;
    }
  else
    return SolverControl::check(step, check_value);
}

/*------------------------ ConsecutiveControl -------------------------------*/


ConsecutiveControl::ConsecutiveControl(const unsigned int n,
                                       const double       tolerance,
                                       const unsigned int n_consecutive_iterations,
                                       const bool m_log_history,
                                       const bool m_log_result)
  :
  SolverControl (n, tolerance, m_log_history, m_log_result),
  n_consecutive_iterations(n_consecutive_iterations),
  n_converged_iterations(0)
{
  AssertThrow(n_consecutive_iterations>0,
              ExcMessage("n_consecutive_iterations should be positive"));
}



ConsecutiveControl::ConsecutiveControl (const SolverControl &c)
  :
  SolverControl(c),
  n_consecutive_iterations(1),
  n_converged_iterations(0)
{}



ConsecutiveControl &
ConsecutiveControl::operator= (const SolverControl &c)
{
  SolverControl::operator=(c);
  n_consecutive_iterations = 1;
  n_converged_iterations = 0;
  return *this;
}



ConsecutiveControl::~ConsecutiveControl()
{}


SolverControl::State
ConsecutiveControl::check (const unsigned int step,
                           const double check_value)
{
  // reset the counter if ConsecutiveControl is being reused
  if (step==0)
    n_converged_iterations = 0;

  SolverControl::State state = SolverControl::check(step, check_value);
  // check if we need to override the success:
  if (state == success)
    {
      n_converged_iterations++;
      if (n_converged_iterations == n_consecutive_iterations)
        {
          return success;
        }
      else
        {
          lcheck = iterate;
          return iterate;
        }
    }
  else
    {
      n_converged_iterations = 0;
      return state;
    }
}

DEAL_II_NAMESPACE_CLOSE
