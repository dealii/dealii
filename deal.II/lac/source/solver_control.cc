// $Id$

#include <base/logstream.h>
#include <lac/solver_control.h>



/*----------------------- SolverControl ---------------------------------*/


SolverControl::SolverControl (const unsigned int maxiter,
			      const double tolerance,
			      const bool log_history) :
		maxsteps(maxiter),
		tol(tolerance),
		lvalue(1.e300),
		lstep(0),
		log_history(log_history)
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
      deallog << "Failure step " << step
	      << " value " << check_value << endl;
      return failure;
    }
  if (check_value <= tol)
    {
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
				   const double red) :
		SolverControl (n, tol),
		reduce(red)
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

