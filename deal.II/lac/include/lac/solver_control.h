//----------------------------  solver_control.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  solver_control.h  ---------------------------
#ifndef __deal2__solver_control_h
#define __deal2__solver_control_h


#include <base/subscriptor.h>


/**
 * Control class for iterative solvers.
 *
 * Used by iterative methods to
 * determine whether the iteration should be continued. To this respect,
 * the virtual function #check()# is called in each iteration
 * with the current iteration
 * step and the value indicating convergence (usually the residual).
 *
 * After the iteration has terminated, the functions #last_value# and
 * #last_step# can be used to obtain information about the final state
 * of the iteration.
 *
 * #check()# can be replaced in derived classes to allow for more
 * sophisticated tests.
 *
 *
 * \section{State}
 * The return states of the check function are of type #State#, which is an
 * enum local to this class. It indicates the state the
 * solver is in.
 *
 * The possible values of State are
 * \begin{itemize}
 * \item #iterate = 0#: continue the iteration.
 * \item #success#: the goal is reached, the iterative method can terminate
 *       successfully.
 * \item #failure#: the iterative method should stop because convergence 
 *       could not be achieved or at least was not achieved within the given
 *       maximal number of iterations.
 * \end{itemize}
 */
class SolverControl : public Subscriptor
{
  public:

				     /**
				      * #Enum# denoting the different
				      * states a solver can be in. See
				      * the general documentation of
				      * this class for more
				      * information.
				      */
    enum State {
      iterate = 0, success, failure
    };
    
				     /**
				      * Constructor. The parameters
				      * #n# and #tol# are the
				      * maximum number of iteration
				      * steps before failure and the
				      * tolerance to determine success
				      * of the iteration.
				      *
				      * #log_history# specifies whether
				      * the history (i.e. the value to
				      * be checked and the number of
				      * the iteration step) shall be
				      * printed to #deallog# stream.
				      * Default is: do not print. Similarly,
				      *  #log_result#
				      * specifies the whether the final result is logged
				      * to #deallog#. Default is yes.
				      */
    SolverControl (const unsigned int n, const double tol,
		   const bool log_history = false,
		   const bool log_result = true);
    
				     /**
				      * Virtual destructor is needed
				      * as there are virtual functions
				      * in this class.
				      */
    virtual ~SolverControl();
    
				     /**
				      * Decide about success or failure
				      * of an iteration.  This function
				      * gets the current iteration step
				      * to determine, whether the
				      * allowed number of steps has
				      * been exceeded and returns
				      * #failure# in this case. If
				      * #check_value# is below the
				      * prescribed tolerance, it
				      * returns #success#. In all
				      * other cases #iterate# is
				      * returned to suggest
				      * continuation of the iterative
				      * procedure.
				      *
				      * #check()# additionally
				      * preserves #step# and
				      * #check_value#. These
				      * values are accessible by
				      * #last_value()# and
				      * #last_step()#.
				      *
				      * Derived classes may overload this
				      * function, e.g. to log the convergence
				      * indicators (#check_value#) or to do
				      * other computations.
				      */
    virtual State check (const unsigned int step, const double check_value);

				     /**
				      * Return the convergence value of last
				      * iteration step for which #check# was
				      * called by the solver.
				      */
    double last_value() const;
    
				     /**
				      * Number of last iteration step.
				      */
    unsigned int last_step() const;

  protected:
				     /**
				      * Maximum number of steps.
				      */
    const unsigned int maxsteps;
    
				     /**
				      * Prescribed tolerance to be achieved.
				      */
    const double       tol;
    
				     /**
				      * Last value of the convergence criterion.
				      */
    double             lvalue;
    
				     /**
				      * Last step.
				      */
    unsigned int       lstep;

				     /**
				      * Log convergence history to #deallog#?
				      */
    const bool         log_history;
				     /**
				      * Log iteration result to #deallog#?
				      * If true, after finishing the iteration, a
				      * statement about failure or success together with
				      * #lstep# and #lvalue# are logged.
				      */
    const bool         log_result;
};


/**
 * Specialization of #SolverControl# which returns #success# if either
 * the specified tolerance is achieved or if the initial residual (or
 * whatever criterion was chosen by the solver class) is reduced by a
 * given factor. This is useful in cases where you don't want to solve
 * exactly, but rather want to gain two digits or if the maximal
 * number of iterations is achived.  For exemple: The maximal number
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
    ReductionControl (const unsigned int maxiter,
		      const double tolerance,
		      const double reduce,
		      const bool log_history = false,
		      const bool log_result = true);

				     /**
				      * Virtual destructor is needed
				      * as there are virtual functions
				      * in this class.
				      */
    virtual ~ReductionControl();
    
				     /**
				      * Decide about success or failure
				      * of an iteration.  This function
				      * calls the one in the base
				      * class, but sets the tolerance
				      * to #reduction * initial value#
				      * upon the first iteration.
				      */
    virtual State check (const unsigned int step,
			 const double check_value);

				     /**
				      * Return the initial convergence
				      * criterion.
				      */
    double initial_value() const;


protected:
				     /**
				      * Desired reduction factor.
				      */
    const double reduce;
    
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


#endif
