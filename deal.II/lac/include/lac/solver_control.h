/*----------------------------   solver_control.h     ---------------------------*/
/*      $Id$                 */
#ifndef __solver_control_H
#define __solver_control_H
/*----------------------------   solver_control.h     ---------------------------*/



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
 */
class SolverControl {
  public:
				     /**
				      * Return states of the check
				      * function, which indicate the state the
				      * solver is in.
				      *
				      * The possible values of State are
				      * <OL>
				      * <LI> #iterate = 0#: continue
				      * the iteration.
				      * <LI> #success#: the goal is reached,
				      * the iterative method can terminate
				      * successfully.
				      * <LI> #failure#!: the iterative
				      * method should stop because
				      * convergence cannot be achieved or at
				      * least was not achieved within the given
				      * maximal number of iterations.
				      * </OL>
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
				      * Default is: do not print.
				      */
    SolverControl (const unsigned int n, const double tol,
		   const bool log_history = false);
    
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
};




/**
 * Specialization  of #SolverControl# which returns #success# if either
 * the specified tolerance is achieved or if the
 * initial residual (or whatever criterion was chosen by the solver
 * class) is reduced by a given factor. This is useful in cases where
 * you don't want to solve exactly, but rather want to gain two digits.
 */
class ReductionControl : public SolverControl {
  public:
				     /**
				      * Constructor.  Provide the
				      * reduction factor additional to
				      * the arguments of the Control
				      * constructor.
				      */
    ReductionControl (const unsigned int maxiter,
		      const double tolerance,
		      const double reduce);
    
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




/*----------------------------   solver_control.h     ---------------------------*/
/* end of #ifndef __solver_control_H */
#endif
/*----------------------------   solver_control.h     ---------------------------*/
