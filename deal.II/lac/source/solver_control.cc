// $Id$

#include <base/logstream.h>
#include <lac/solver_control.h>


/*----------------------- SolverControl ---------------------------------*/


SolverControl::SolverControl (const unsigned int maxiter,
			      const double tolerance,
			      const bool log_history,
			      const bool log_result)
		:
		maxsteps(maxiter),
		tol(tolerance),
		lvalue(1.e300),
		lstep(0),
		log_history(log_history),
		log_result(log_result)
{};


SolverControl::~SolverControl() 
{};


SolverControl::State
SolverControl::check (const unsigned int step,
		      const double check_value)
{
  if (log_history)
    deallog << "Check " << step << "\t" << check_value << endl;
  
  lstep  = step;
  lvalue = check_value;
  if (step>=maxsteps)
    {
      if (log_result)
	deallog << "Failure step " << step
		<< " value " << check_value << endl;
      return failure;
    }
  if (check_value <= tol)
    {
      if (log_result)
	deallog << "Convergence step " << step
		<< " value " << check_value << endl;      
      return success;
    }
  return iterate;
};



double
SolverControl::last_value() const
{
  return lvalue;
};



unsigned int
SolverControl::last_step() const
{
  return lstep;
};




/*----------------------- ReductionControl ---------------------------------*/



ReductionControl::ReductionControl(const unsigned int n,
				   const double tol,
				   const double red,
				   const bool log_history,
				   const bool log_result)
		:
		SolverControl (n, tol, log_history, log_result),
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
    return success;
  else
    return SolverControl::check(step, check_value);
};

