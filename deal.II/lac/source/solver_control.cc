//----------------------------  solver_control.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002 by the deal.II authors
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

#ifdef HAVE_STD_STRINGSTREAM
#  include <sstream>
#else
#  include <strstream>
#endif

/*----------------------- SolverControl ---------------------------------*/


SolverControl::NoConvergence::NoConvergence (const unsigned int last_step,
					     const double       last_residual)
		:
		last_step (last_step),
		last_residual (last_residual)
{};


const char *
SolverControl::NoConvergence::what () const throw ()
{
  				   // have a place where to store the
				   // description of the exception as a char *
				   //
				   // this thing obviously is not multi-threading
				   // safe, but we don't care about that for now
				   //
				   // we need to make this object static, since
				   // we want to return the data stored in it
				   // and therefore need a liftime which is
				   // longer than the execution time of this
				   // function
  static std::string description;
				   // convert the messages printed by the
				   // exceptions into a std::string
#ifdef HAVE_STD_STRINGSTREAM
  std::ostringstream out;
#else
  std::ostrstream out;
#endif

  out << "Iterative method reported convergence failure in step "
      << last_step << " with residual " << last_residual;
#ifndef HAVE_STD_STRINGSTREAM
  out << std::ends;
#endif

  description = out.str();
  return description.c_str();
}



SolverControl::SolverControl (const unsigned int maxiter,
			      const double tolerance,
			      const bool m_log_history,
			      const bool m_log_result)
		:
		maxsteps(maxiter),
		tol(tolerance),
		lvalue(1.e300),
		lstep(0),
		check_failure(false),
		relative_failure_residual(0),
		failure_residual(0),
		m_log_history(m_log_history),
		m_log_frequency(1),
		m_log_result(m_log_result)
{};



SolverControl::~SolverControl() 
{};



SolverControl::State
SolverControl::check (const unsigned int step,
		      const double check_value)
{
  if (m_log_history && ((step % m_log_frequency) == 0))
    deallog << "Check " << step << "\t" << check_value << std::endl;
  
  lstep  = step;
  lvalue = check_value;

  if (step==0)
    {
      if (check_failure)
	failure_residual=relative_failure_residual*check_value;
      
      if (m_log_result)
	deallog << "Starting value " << check_value << std::endl;
    }

  if (check_value <= tol)
    {
      if (m_log_result)
	deallog << "Convergence step " << step
		<< " value " << check_value << std::endl;
      lcheck = success;
      return success;
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
      if (m_log_result)
	deallog << "Failure step " << step
		<< " value " << check_value << std::endl;
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
  unsigned int old = m_log_frequency;
  m_log_frequency = f;
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
				   const bool m_log_history,
				   const bool m_log_result)
		:
		SolverControl (n, tol, m_log_history, m_log_result),
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
				   // if this is the first time we
				   // come here, then store the
				   // residual for later comparisons
  if (step==0)
    {
      initial_val = check_value;
      reduced_tol = check_value * reduce;
    };

				   // check whether desired reduction
				   // has been achieved. also check
				   // for equality in case initial
				   // residual already was zero
  if (check_value <= reduced_tol)
    {
      if (m_log_result)
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
