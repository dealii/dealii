//----------------------------  solver_control.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  solver_control.cc  ---------------------------


#include <base/logstream.h>
#include <base/parameter_handler.h>
#include <lac/solver_control.h>

#include <cmath>

/*----------------------- SolverControl ---------------------------------*/


SolverControl::SolverControl (const unsigned int maxiter,
			      const double tolerance,
			      const bool _log_history,
			      const bool _log_result)
		:
		maxsteps(maxiter),
		tol(tolerance),
		lvalue(1.e300),
		lstep(0),
		check_failure(false),
		relative_failure_residual(0),
		failure_residual(0),
		_log_history(_log_history),
		_log_frequency(1),
		_log_result(_log_result)
{};



SolverControl::~SolverControl() 
{};



SolverControl::State
SolverControl::check (const unsigned int step,
		      const double check_value)
{
  if (_log_history && ((step % _log_frequency) == 0))
    deallog << "Check " << step << "\t" << check_value << std::endl;
  
  lstep  = step;
  lvalue = check_value;

  if (step==0)
    {
      if (check_failure)
	failure_residual=relative_failure_residual*check_value;
      
      if (_log_result)
	deallog << "Starting value " << check_value << std::endl;
    }

  
  if ((step >= maxsteps) ||
#ifdef HAVE_ISNAN
      isnan(check_value) ||
#else
#  if HAVE_UNDERSCORE_ISNAN
				       // on Microsoft Windows, the
				       // function is called _isnan
      _isnan(check_value) ||
#  endif
#endif
      (check_failure && (check_value > failure_residual))
  )
    {
      if (_log_result)
	deallog << "Failure step " << step
		<< " value " << check_value << std::endl;
      lcheck = failure;
      return failure;
    }

  if (check_value <= tol)
    {
      if (_log_result)
	deallog << "Convergence step " << step
		<< " value " << check_value << std::endl;
      lcheck = success;
      return success;
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
SolverControl::last_value() const
{
  return lvalue;
}


unsigned int
SolverControl::last_step() const
{
  return lstep;
}


unsigned int
SolverControl::log_frequency (unsigned int f)
{
  if (f==0)
    f = 1;
  unsigned int old = _log_frequency;
  _log_frequency = f;
  return old;
}


void
SolverControl::declare_parameters (ParameterHandler& param)
{
  param.declare_entry ("Max steps", "100", Patterns::Integer());
  param.declare_entry ("Tolerance", "1.e-10", Patterns::Double());
  param.declare_entry ("Log history", "false", Patterns::Bool());
  param.declare_entry ("Log frequency", "1", Patterns::Integer());
  param.declare_entry ("Log result", "true", Patterns::Bool());
}


void SolverControl::parse_parameters (ParameterHandler& param)
{
  set_max_steps (param.get_integer("Max steps"));
  set_tolerance (param.get_double("Tolerance"));
  log_history (param.get_bool("Log history"));
  log_result (param.get_bool("Log result"));
  log_frequency (param.get_integer("Log frequency"));
}

/*----------------------- ReductionControl ---------------------------------*/


ReductionControl::ReductionControl(const unsigned int n,
				   const double tol,
				   const double red,
				   const bool _log_history,
				   const bool _log_result)
		:
		SolverControl (n, tol, _log_history, _log_result),
		reduce(red)
{};


ReductionControl::~ReductionControl()
{};


double
ReductionControl::initial_value() const {
  return initial_val;
};


SolverControl::State
ReductionControl::check (const unsigned int step,
			 const double check_value)
{
  if (step==0)
    {
      initial_val = check_value;
      reduced_tol = check_value * reduce;
    };

  if (check_value < reduced_tol)
    {
      if (_log_result)
	deallog << "Convergence step " << step
		<< " value " << check_value << std::endl;
      lcheck = success;
      return success;
    }
  else
    return SolverControl::check(step, check_value);
};



void
ReductionControl::declare_parameters (ParameterHandler& param)
{
  SolverControl::declare_parameters (param);
  param.declare_entry("Reduction", "1.e-2", Patterns::Double());
}


void
ReductionControl::parse_parameters (ParameterHandler& param)
{
  SolverControl::parse_parameters (param);
  set_reduction (param.get_double("Reduction"));
}

