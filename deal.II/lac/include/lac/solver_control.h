//----------------------------  solver_control.h  ---------------------------
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
//----------------------------  solver_control.h  ---------------------------
#ifndef __deal2__solver_control_h
#define __deal2__solver_control_h


#include <base/config.h>
#include <base/subscriptor.h>

class ParameterHandler;


/**
 * Control class for iterative solvers.
 *
 * Used by iterative methods to determine whether the iteration should
 * be continued. To this respect, the virtual function @p{check()} is
 * called in each iteration with the current iteration step and the
 * value indicating convergence (usually the residual).
 *
 * After the iteration has terminated, the functions @p{last_value} and
 * @p{last_step} can be used to obtain information about the final state
 * of the iteration.
 *
 * @p{check()} can be replaced in derived classes to allow for more
 * sophisticated tests.
 *
 *
 * @sect2{State}
 * The return states of the check function are of type @p{State},
 * which is an enum local to this class. It indicates the state the
 * solver is in.
 *
 * The possible values of State are
 * @begin{itemize}
 * @item @p{iterate = 0}: continue the iteration.
 * @item @p{success}: the goal is reached, the iterative method can terminate
 *       successfully.
 * @item @p{failure}: the iterative method should stop because convergence 
 *       could not be achieved or at least was not achieved within the given
 *       maximal number of iterations.
 * @end{itemize}
 */
class SolverControl : public Subscriptor
{
  public:

				     /**
				      * @p{Enum} denoting the different
				      * states a solver can be in. See
				      * the general documentation of
				      * this class for more
				      * information.
				      */
    enum State {
      iterate = 0, success, failure
    };

				     /**
				      * Class to be thrown upon
				      * failing convergence of an
				      * iterative solver, when either
				      * the number of iterations
				      * exceeds the limit or the
				      * residual fails to reach the
				      * desired limit, e.g. in the
				      * case of a break-down.
				      *
				      * The residual in the last
				      * iteration, as well as the
				      * iteration number of the last
				      * step are stored in this object
				      * and can be recovered upon
				      * catching an exception of this
				      * class.
				      */
    class NoConvergence : public std::exception 
    {
      public:
					 /**
					  * Constructor.
					  */
	NoConvergence (const unsigned int last_step,
		       const double       last_residual);

					 /**
					  * Standardized output for
					  * catch handlers.
					  */
	virtual const char * what () const throw ();

					 /**
					  * Iteration number of the
					  * last step.
					  */
	const unsigned int last_step;

					 /**
					  * Residual in the last step.
					  */
	const double       last_residual;
    };
    
    
				     /**
				      * Constructor. The parameters
				      * @p{n} and @p{tol} are the
				      * maximum number of iteration
				      * steps before failure and the
				      * tolerance to determine success
				      * of the iteration.
				      *
				      * @p{log_history} specifies
				      * whether the history (i.e. the
				      * value to be checked and the
				      * number of the iteration step)
				      * shall be printed to
				      * @p{deallog} stream.  Default
				      * is: do not print. Similarly,
				      * @p{log_result} specifies the
				      * whether the final result is
				      * logged to @p{deallog}. Default
				      * is yes.
				      */
    SolverControl (const unsigned int n           = 100,
		   const double       tol         = 1.e-10,
		   const bool         log_history = false,
		   const bool         log_result  = true);
    
				     /**
				      * Virtual destructor is needed
				      * as there are virtual functions
				      * in this class.
				      */
    virtual ~SolverControl();
    
				     /**
				      * Interface to parameter file.
				      */
    static void declare_parameters (ParameterHandler& param);

				     /**
				      * Read parameters from file.
				      */
    void parse_parameters (ParameterHandler& param);

				     /**
				      * Decide about success or failure
				      * of an iteration.  This function
				      * gets the current iteration step
				      * to determine, whether the
				      * allowed number of steps has
				      * been exceeded and returns
				      * @p{failure} in this case. If
				      * @p{check_value} is below the
				      * prescribed tolerance, it
				      * returns @p{success}. In all
				      * other cases @p{iterate} is
				      * returned to suggest
				      * continuation of the iterative
				      * procedure.
				      *
				      * The iteration is also aborted
				      * if the residual becomes a
				      * denormalized value
				      * (@p{NaN}). Note, however, that
				      * this check is only performed
				      * if the @p{isnan} function is
				      * provided by the operating
				      * system, which is not always
				      * true. The @p{configure}
				      * scripts checks for this and
				      * sets the flag @p{HAVE_ISNAN}
				      * in the file
				      * @p{Make.global_options} if
				      * this function was found.
				      *
				      * @p{check()} additionally
				      * preserves @p{step} and
				      * @p{check_value}. These
				      * values are accessible by
				      * @p{last_value()} and
				      * @p{last_step()}.
				      *
				      * Derived classes may overload
				      * this function, e.g. to log the
				      * convergence indicators
				      * (@p{check_value}) or to do
				      * other computations.
				      */
    virtual State check (const unsigned int step,
			 const double   check_value);

				     /**
				      * Return the result of the last check operation.
				      */
    State last_check() const;
    
				     /**
				      * Return the convergence value of last
				      * iteration step for which @p{check} was
				      * called by the solver.
				      */
    double last_value() const;
    
				     /**
				      * Number of last iteration step.
				      */
    unsigned int last_step() const;

				     /**
				      * Maximum number of steps.
				      */
    unsigned int max_steps () const;
    
				     /**
				      * Change maximum number of steps.
				      */
    unsigned int set_max_steps (const unsigned int);

				     /**
				      * Enables the failure
				      * check. Solving is stopped with
				      * @p{ReturnState} @p{failure} if
				      * @p{residual>failure_residual} with
				      * @p{failure_residual:=rel_failure_residual*first_residual}.
				      */
    void set_failure_criterion (const double rel_failure_residual);

				     /**
				      * Disables failure check and
				      * resets
				      * @p{relative_failure_residual}
				      * and @p{failure_residual} to
				      * zero.
				      */
    void clear_failure_criterion ();
    
				     /**
				      * Tolerance.
				      */
    double tolerance () const;

				     /**
				      * Change tolerance.
				      */
    double set_tolerance (const double);

				     /**
				      * Log each iteration step. Use
				      * @p{log_frequency} for skipping
				      * steps.
				      */
    void log_history (const bool);

				     /**
				      * Returns the log_history flag.
				      */
    bool log_history () const;
    
				     /**
				      * Set logging frequency.
				      */
    unsigned int log_frequency (unsigned int);

				     /**
				      * Log start and end step.
				      */
    void log_result (const bool);
    
  protected:
				     /**
				      * Maximum number of steps.
				      */
    unsigned int maxsteps;
    
				     /**
				      * Prescribed tolerance to be achieved.
				      */
    double       tol;

				     /**
				      * Result of last check operation.
				      */
    State        lcheck;
    
				     /**
				      * Last value of the convergence criterion.
				      */
    double       lvalue;
    
				     /**
				      * Last step.
				      */
    unsigned int lstep;

				     /**
				      * Is set to @p{true} by
				      * @p{set_failure_criterion} and
				      * enables failure checking.
				      */
    bool         check_failure;

				     /*
				      * Stores the
				      * @p{rel_failure_residual} set by
				      * @p{set_failure_criterion}
				      */
    double       relative_failure_residual;

    
				     /**
				      * @p{failure_residual} equals the
				      * first residual multiplied by
				      * @p{relative_crit} set by
				      * @p{set_failure_criterion} (see there).
				      *
				      * Until the first residual is
				      * known it is 0.
				      */
    double       failure_residual;
    
				     /**
				      * Log convergence history to
				      * @p{deallog}.
				      */
    bool         m_log_history;
				     /**
				      * Log only every nth step.
				      */
    unsigned int m_log_frequency;
    
				     /**
				      * Log iteration result to
				      * @p{deallog}.  If true, after
				      * finishing the iteration, a
				      * statement about failure or
				      * success together with @p{lstep}
				      * and @p{lvalue} are logged.
				      */
    bool         m_log_result;
};


/**
 * Specialization of @p{SolverControl} which returns @p{success} if either
 * the specified tolerance is achieved or if the initial residual (or
 * whatever criterion was chosen by the solver class) is reduced by a
 * given factor. This is useful in cases where you don't want to solve
 * exactly, but rather want to gain two digits or if the maximal
 * number of iterations is achieved.  For example: The maximal number
 * of iterations is 20, the reduction factor is 1% und the tolerance
 * is 0.1%. The initial residual is 2.5. The process will break if 20
 * iteration are comleted or the new residual is less then 2.5*1% or
 * if it is less then 0.1%. 
 */
class ReductionControl : public SolverControl
{
  public:
				     /**
				      * Constructor.  Provide the
				      * reduction factor additional to
				      * the arguments of the Control
				      * constructor.
				      */
    ReductionControl (const unsigned int maxiter = 100,
		      const double   tolerance   = 1.e-10,
		      const double   reduce      = 1.e-2,
		      const bool     log_history = false,
		      const bool     log_result  = true);

				     /**
				      * Virtual destructor is needed
				      * as there are virtual functions
				      * in this class.
				      */
    virtual ~ReductionControl();
    
				     /**
				      * Interface to parameter file.
				      */
    static void declare_parameters (ParameterHandler& param);

				     /**
				      * Read parameters from file.
				      */
    void parse_parameters (ParameterHandler& param);
    
				     /**
				      * Decide about success or failure
				      * of an iteration.  This function
				      * calls the one in the base
				      * class, but sets the tolerance
				      * to @p{reduction * initial value}
				      * upon the first iteration.
				      */
    virtual State check (const unsigned int step,
			 const double   check_value);

				     /**
				      * Return the initial convergence
				      * criterion.
				      */
    double initial_value() const;

				     /**
				      * Reduction factor.
				      */
    double reduction () const;

				     /**
				      * Change reduction factor.
				      */
    double set_reduction (const double);

  protected:
				     /**
				      * Desired reduction factor.
				      */
    double reduce;
    
				     /**
				      * Initial value.
				      */
    double initial_val;
    
				     /**
				      * Reduced tolerance. Stop iterations
				      * if either this value is achieved
				      * or if the base class indicates
				      * success.
				      */
    double reduced_tol;
};

//------------------------------------------------------------//


inline unsigned int
SolverControl::max_steps () const
{
  return maxsteps;
}



inline unsigned int
SolverControl::set_max_steps (const unsigned int newval)
{
  unsigned int old = maxsteps;
  maxsteps = newval;
  return old;
}



inline void
SolverControl::set_failure_criterion (const double rel_failure_residual)
{
  relative_failure_residual=rel_failure_residual;
  check_failure=true;
}



inline void
SolverControl::clear_failure_criterion ()
{
  relative_failure_residual=0;
  failure_residual=0;
  check_failure=false;
}



inline double
SolverControl::tolerance () const
{
  return tol;
}



inline double
SolverControl::set_tolerance (const double t)
{
  double old = tol;
  tol = t;
  return old;
}


inline void
SolverControl::log_history (const bool newval)
{
  m_log_history = newval;
}



inline bool
SolverControl::log_history () const
{
  return m_log_history;
}



inline void
SolverControl::log_result (const bool newval)
{
  m_log_result = newval;
}



inline double
ReductionControl::reduction () const
{
  return reduce;
}



inline double
ReductionControl::set_reduction (const double t)
{
  double old = reduce;
  reduce = t;
  return old;
}

#endif
