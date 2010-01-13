//---------------------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2006, 2007, 2010 by Guido Kanschat
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#include <numerics/time_step_control.h>
#include <base/parameter_handler.h>

DEAL_II_NAMESPACE_OPEN

using namespace Algorithms;

TimestepControl::TimestepControl (double start,
				  double final,
				  double tolerance,
				  double start_step,
				  double print_step,
				  double max_step)
  : start_val(start),
    final_val(final),
    tolerance_val(tolerance),
    strategy_val(uniform),
    start_step_val(start_step),
    max_step_val(max_step),
    min_step_val(0),
    current_step_val(start_step),
    step_val(start_step),
    print_step(print_step)
{
  now_val = start_val;
  strcpy(format, "T.%06.3f");
}



void TimestepControl::declare_parameters (ParameterHandler& param)
{
  param.declare_entry ("Start", "0.", Patterns::Double());
  param.declare_entry ("Final", "1.", Patterns::Double());
  param.declare_entry ("First step", "1.e-2", Patterns::Double());
  param.declare_entry ("Max step", "1.", Patterns::Double());
  param.declare_entry ("Tolerance", "1.e-2", Patterns::Double());
  param.declare_entry ("Print step", "-1.", Patterns::Double());
  param.declare_entry ("Strategy", "uniform",
		       Patterns::Selection("uniform|doubling"));
}




void TimestepControl::parse_parameters (ParameterHandler& param)
{
  start (param.get_double ("Start"));
  start_step (param.get_double ("First step"));
  max_step (param.get_double ("Max step"));
  final (param.get_double ("Final"));
  tolerance (param.get_double ("Tolerance"));
  print_step = param.get_double ("Print step");
  const std::string strat = param.get("Strategy");
  if (strat == std::string("uniform"))
    strategy_val = uniform;
  else if (strat == std::string("doubling"))
    strategy_val = doubling;
}




bool
TimestepControl::advance ()
{
  bool changed = false;
  double s = step_val;

				   // Do time step control, but not in
				   // first step.
  if (now_val != start())
    {
      if (strategy_val == doubling && 2*s <= tolerance_val)
	s *= 2;
      if (s > max_step_val)
	s = max_step_val;
    }

				   // Try incrementing time by s
  double h = now_val + s;
  changed = s != step_val;
  
  step_val = s;
  current_step_val = s;
				   // If we just missed the final
				   // time, increase the step size a
				   // bit. This way, we avoid a very
				   // small final step. If the step
				   // shot over the final time, adjust
				   // it so we hit the final time
				   // exactly.
  double s1 = .01*s;
  if (h > final_val-s1)
    {
      current_step_val = final_val - now_val;
      h = final_val;
      changed = true;
    }
  
  now_val = h;
  return changed;
}


bool TimestepControl::print ()
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

