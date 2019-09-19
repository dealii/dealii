// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2019 by the deal.II authors
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


#include <deal.II/algorithms/timestep_control.h>

#include <deal.II/base/parameter_handler.h>

DEAL_II_NAMESPACE_OPEN

using namespace Algorithms;

TimestepControl::TimestepControl(double start,
                                 double final,
                                 double tolerance,
                                 double start_step,
                                 double print_step,
                                 double max_step)
  : start_val(start)
  , final_val(final)
  , tolerance_val(tolerance)
  , start_step_val(start_step)
  , max_step_val(max_step)
  , min_step_val(0)
  , current_step_val(start_step)
  , step_val(start_step)
  , print_step(print_step)
  , next_print_val(print_step > 0. ? start_val + print_step : start_val - 1.)
{
  now_val = start_val;

  // avoid compiler warning
  (void)min_step_val;
}



void
TimestepControl::declare_parameters(ParameterHandler &param)
{
  param.declare_entry("Start", "0.", Patterns::Double());
  param.declare_entry("Final", "1.", Patterns::Double());
  param.declare_entry("First step", "1.e-2", Patterns::Double(0.));
  param.declare_entry("Max step", "1.", Patterns::Double(0.));
  param.declare_entry("Tolerance", "1.e-2", Patterns::Double(0.));
  param.declare_entry("Print step", "-1.", Patterns::Double());
}



void
TimestepControl::parse_parameters(ParameterHandler &param)
{
  start(param.get_double("Start"));
  start_step(param.get_double("First step"));
  max_step(param.get_double("Max step"));
  final(param.get_double("Final"));
  tolerance(param.get_double("Tolerance"));
  print_step = param.get_double("Print step");
}



bool
TimestepControl::advance()
{
  bool changed = false;

  // Try incrementing time by s
  double now_trial = now_val + step_val;
  current_step_val = step_val;

  // If we just missed the final time, increase the step size a bit. This way,
  // we avoid a very small final step. If the step shot over the final time,
  // adjust it so we hit the final time exactly.
  double s1 = .01 * step_val;
  if (now_trial > final_val - s1)
    {
      current_step_val = final_val - now_val;
      now_trial        = final_val;
      changed          = true;
    }

  now_val = now_trial;
  return changed;
}


bool
TimestepControl::print()
{
  if (print_step == 0.)
    return false;
  if (print_step < 0.)
    return true;

  bool result = (now_val >= next_print_val);

  if (result)
    {
      next_print_val += print_step;
      if (next_print_val > final_val)
        next_print_val = final_val;
    }
  return result;
}

DEAL_II_NAMESPACE_CLOSE
