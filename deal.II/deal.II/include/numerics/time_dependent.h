/*----------------------------   time-dependent.h     ---------------------------*/
/*      $Id$                 */
#ifndef __time_dependent_H
#define __time_dependent_H
/*----------------------------   time-dependent.h     ---------------------------*/


#include <base/exceptions.h>
#include <base/subscriptor.h>
#include <basic/forward-declarations.h>

#include <vector>



class TimeDependent 
{
  public:
    struct TimeSteppingData
    {
					 /**
					  * Constructor; see the different
					  * fields for a description of the
					  * meaning of the parameters.
					  */
	TimeSteppingData (const unsigned int look_ahead,
			  const unsigned int look_back);
	
					 /**
					  * This denotes the number of timesteps
					  * the timestepping algorithm needs to
					  * look ahead. Usually, this number
					  * will be zero, since algorithms
					  * looking ahead can't act as
					  * timestepping schemes since they
					  * can't compute their data from knowledge
					  * of the past only and are therefore
					  * global in time.
					  *
					  * However, it may be necessary to look
					  * ahead in other circumstances, when
					  * not wanting to access the data of the
					  * next time step(s), but for example
					  * to know the next grid, the solution
					  * of a dual problem on the next
					  * time level, etc.
					  *
					  * Note that for a dual problem walking
					  * back in time, "looking ahead" means
					  * looking towards smaller time values.
					  *
					  * The value of this number determines,
					  * how many time steps ahead the
					  * time step manager start to call
					  * the #wake_up# function for each
					  * time step.
					  */
	const unsigned int look_ahead;
	
					 /**
					  * This is the opposite variable to the
					  * above one. It denotes the number of
					  * time steps behind the present one
					  * for which we need to keep all data
					  * in order to do the computations on
					  * the present time level.
					  *
					  * For one step schemes (e.g. the
					  * Euler schemes, or the Crank-Nicolson
					  * scheme), this value will be one.
					  *
					  * The value of this number
					  * determines, how many time
					  * steps after having done
					  * computations on a tim level
					  * the time step manager will
					  * call the #sleep# function for
					  * each time step.
					  */
	const unsigned int look_back;
    };
    

				     /**
				      * Constructor.
				      */
    TimeDependent (const TimeSteppingData &data_primal,
		   const TimeSteppingData &data_dual,
		   const TimeSteppingData &data_postprocess);
    
    
				     /**
				      * Destructor. This will delete the
				      * objects pointed to by the pointers
				      * given to the #insert_*# and
				      * #add_timestep# functions, i.e.
				      * it will delete the objects doing
				      * the computations on each time step.
				      */
    virtual ~TimeDependent ();

				     /**
				      * Add a timestep at any position. The
				      * position may be zero (at the start)
				      * through #N# (at the end), where
				      * #N# is the number of timesteps
				      * stored in this object previously.
				      *
				      * Note that by giving an object
				      * to this function, the
				      * #TimeDependent# object assumes
				      * ownership of the object; it will
				      * therefore also take care of
				      * deletion of the objects its manages.
				      * This mechanism usually will result
				      * in a set-up loop like this
				      * \begin{verbatim}
				      * for (i=0; i<N; ++i)
				      *   manager.add_timestep(new MyTimeStep());
				      * \end{verbatim}
				      *
				      * There is another function,
				      * #add_timestep#, which inserts a
				      * time step at the end of the list.
				      */
    void insert_timestep (TimeStepBase      *new_timestep,
			  const unsigned int position);

				     /**
				      * Just like #insert_timestep#, but
				      * insert at the end.
				      */
    void add_timestep (TimeStepBase *new_timestep);

				     /**
				      * Delete a timestep. This is only
				      * necessary to call, if you want to
				      * delete it between two sweeps; at
				      * the end of the lifetime of this object,
				      * care is taken automatically of
				      * deletion of the time step objects.
				      * Deletion of the object by the
				      * destructor is done through this
				      * function also.
				      */
    void delete_timestep (const unsigned int position);
    
				     /**
				      * Solve the primal problem; uses the
				      * functions #init_for_primal_problem#
				      * and #solve_primal_problem# of the
				      * #TimeStepBase# class through the
				      * #do_loop# function of this class.
				      *
				      * Look ahead and look back are
				      * determined by the #timestepping_data_primal#
				      * object given to the constructor.
				      */
    void solve_primal_problem ();

				     /**
				      * Solve the dual problem; uses the
				      * functions #init_for_dual_problem#
				      * and #solve_dual_problem# of the
				      * #TimeStepBase# class through the
				      * #do_loop# function of this class.
				      *
				      * Look ahead and look back are
				      * determined by the #timestepping_data_dual#
				      * object given to the constructor.
				      */
    void solve_dual_problem ();

				     /**
				      * Do a postprocessing round; uses the
				      * functions #init_for_postprocessing#
				      * and #postprocess# of the
				      * #TimeStepBase# class through the
				      * #do_loop# function of this class.
				      *
				      * Look ahead and look back are
				      * determined by the #timestepping_data_postprocess#
				      * object given to the constructor.
				      */
    void postprocess ();
    
				     /**
				      * Do a loop over all timesteps, call the
				      * #init_function# at the beginning and
				      * the #loop_function# of each time step.
				      * The #timestepping_data# determine how
				      * many timesteps in front and behind
				      * the present one the #wake_up# and
				      * #sleep# functions are called.
				      *
				      * To see how this function work, note that
				      * the function #solve_primal_problem# only
				      * consists of a call to
				      * #  do_loop (&TimeStepBase::init_for_primal_problem,
				      *	   &TimeStepBase::solve_primal_problem,
				      *	   timestepping_data_primal);#.
				      *
				      * Note also, that the given class from which
				      * the two functions are taken needs not
				      * necessarily be a #TimeStepBase#, but it
				      * should be a derived class, that is
				      * #static_cast#able from a #TimeStepBase#.
				      * The function may be a virtual function
				      * (even a pure one) of that class, which
				      * should help if the actual class where it
				      * is implemented is one which is derived
				      * through virtual base classes and thus
				      * unreachable by #static_cast# from the
				      * #TimeStepBase# class.
				      */
    template <class TimeStepClass>
    void do_loop (void (TimeStepClass::*init_function) (),
		  void (TimeStepClass::*loop_function) (),
		  const TimeSteppingData &timestepping_data);
    
		  
				     /**
				      * Initialize the objects for the next
				      * sweep. This function specifically does
				      * the following: assign each time
				      * level the number it presently has
				      * within the array (which may change,
				      * if time levels are inserted or
				      * deleted) and transmit the number of
				      * the present sweep to these objects.
				      *
				      * It also calls the #init_for_sweep#
				      * function of each time step object,
				      * after the numbers above are set.
				      *
				      * This function is virtual, so you
				      * may overload it. You should, however
				      * not forget to call this function as
				      * well from your overwritten version,
				      * at best at the beginning of your
				      * function since this is some kind of
				      * "constructor-like" function, which
				      * should be called bottom-up.
				      */
    virtual void start_sweep (const unsigned int sweep_no);

				     /**
				      * Exception.
				      */
    DeclException2 (ExcInvalidPosition,
		    int, int,
		    << "Can\'t insert time step at position " << arg1
                    << " since there only " << arg2 << " positions at all.");
    
  protected:
				     /**
				      * Vector holding pointers to the time
				      * level objects. This is the main data
				      * this object operates on. Note that
				      * this object takes possession of the
				      * objects pointed to by the pointers
				      * in this collection.
				      */
    vector<TimeStepBase*> timesteps;

				     /**
				      * Number of the present sweep. This is
				      * reset by the #start_sweep# function
				      * called at the outset of each sweep.
				      */
    unsigned int sweep_no;

				     /**
				      * Some flags telling the
				      * #solve_primal_problem# function what to
				      * do. See the documentation of this struct
				      * for more information.
				      */
    TimeSteppingData timestepping_data_primal;

				     /**
				      * Some flags telling the
				      * #solve_dual_problem# function what to
				      * do. See the documentation of this struct
				      * for more information.
				      */
    TimeSteppingData timestepping_data_dual;

				     /**
				      * Some flags telling the
				      * #postprocess# function what to
				      * do. See the documentation of this struct
				      * for more information.
				      */
    TimeSteppingData timestepping_data_postprocess;  
    
};




template <class TimeStepClass>
void TimeDependent::do_loop (void (TimeStepClass::*init_function) (),
			     void (TimeStepClass::*loop_function) (),
			     const TimeSteppingData &timestepping_data)
{
  const unsigned int n_timesteps = timesteps.size();

				   // initialize the time steps for
				   // a round of primal problems
  for (unsigned int step=0; step<n_timesteps; ++step)
    (static_cast<TimeStepClass*>(timesteps[step])->*init_function) ();

				   // wake up the first few time levels
  for (int step=-timestepping_data.look_ahead; step<0; ++step)
    for (int look_ahead=0;
	 look_ahead<=static_cast<int>(timestepping_data.look_ahead); ++look_ahead)
      if (step+look_ahead >= 0)
	timesteps[step+look_ahead]->wake_up(look_ahead);
  
  for (unsigned int step=0; step<n_timesteps; ++step)
    {
				       // first thing: wake up the
				       // timesteps ahead as necessary
      for (unsigned int look_ahead=0;
	   look_ahead<=timestepping_data.look_ahead; ++look_ahead)
	if (step+look_ahead < n_timesteps)
	  timesteps[step+look_ahead]->wake_up(look_ahead);
      
				       // actually do the work
      (static_cast<TimeStepClass*>(timesteps[step])->*loop_function) ();
      
				       // let the timesteps behind sleep
      for (unsigned int look_back=0;
	   look_back<=timestepping_data.look_back; ++look_back)
	if (step>=look_back)
	  timesteps[step-look_back]->sleep(look_back);
    };

				   // make the last few timesteps sleep
  for (int step=n_timesteps;
       step<static_cast<int>(n_timesteps+timestepping_data.look_back); ++step)
    for (int look_back=0;
	 look_back<=static_cast<int>(timestepping_data.look_back); ++look_back)
      if ((step-look_back>=0) && (step-look_back<static_cast<int>(n_timesteps)))
	timesteps[step-look_back]->sleep(look_back);
};









class TimeStepBase : public Subscriptor
{
  public:
				     /**
				      * Enum denoting the type of problem
				      * which will have to be solved next.
				      */
    enum NextAction {
	  primal_problem,
	  dual_problem,
	  postprocess
    };
    
				     /**
				      * Constructor. Does nothing here apart
				      * from setting the time.
				      */
    TimeStepBase (const double time);

				     /**
				      * Destructor. At present, this does
				      * nothing.
				      */
    virtual ~TimeStepBase ();

				     /**
				      * Reconstruct all the data that is
				      * needed for this time level to work.
				      * This function serves to reget all
				      * the variables and data structures
				      * to work again after they have been
				      * send to sleep some time before, or
				      * at the first time we visit this time
				      * level. In particular, it is used
				      * to reconstruct the triangulation,
				      * degree of freedom handlers, to reload
				      * data vectors in case they have been
				      * stored to disk, etc.
				      *
				      * The actual implementation of
				      * this function does nothing.
				      *
				      * Since this is an important task, you
				      * should call this function from your
				      * own function, should you choose to
				      * overload it in your own class (which
				      * likely is the case), preferably at
				      * the beginning so that your function
				      * can take effect of the triangulation
				      * already existing.
				      */
    virtual void wake_up (const unsigned);

				     /**
				      * This is the opposite function
				      * to #wake_up#. It is used to
				      * delete data or save it to disk
				      * after they are no more needed
				      * for the present sweep. Typical
				      * kinds of data for this are
				      * data vectors, degree of
				      * freedom handlers,
				      * triangulation objects,
				      * etc. which occupy large
				      * amounts of memory and may
				      * therefore be externalized.
				      *
				      * By default, this function does
				      * nothing.
				      */
    virtual void sleep (const unsigned);

				     /**
				      * This function is called each time
				      * before a new sweep is started. You
				      * may want to set up some fields needed
				      * in the course of the computations,
				      * and so on. You should take good care,
				      * however, not to install large objects,
				      * which should be deferred until the
				      * #wake_up# function is called.
				      *
				      * A typical action of this function
				      * would be sorting out names of
				      * temporary files needed in the
				      * process of solving, etc.
				      *
				      * At the time this function is called,
				      * the values of #timestep_no#, #sweep_no#
				      * and the pointer to the previous and
				      * next time step object already have
				      * their correct value.
				      *
				      * The default implementation of this
				      * function does nothing.
				      */
    virtual void init_for_sweep ();

				     /**
				      * Before the primal problem is
				      * solved on each time level, this
				      * function is called (i.e. before the
				      * solution takes place on the first
				      * time level). By default, this function
				      * sets the #next_action# variable of
				      * this class. You may overload this
				      * function, but you should call this
				      * function within your own one.
				      */
    virtual void init_for_primal_problem ();

				     /**
				      * Same as above, but called before
				      * a round of dual problem solves.
				      */
    virtual void init_for_dual_problem ();

    				     /**
				      * Same as above, but called before
				      * a round of postprocessing steps.
				      */
    virtual void init_for_postprocessing ();

				     /**
				      * This function is called by the
				      * manager object when solving the
				      * primal problem on this time level
				      * is needed. It is called after
				      * the #wake_up# function was
				      * called and before the #sleep#
				      * function will be called. There
				      * is no default implementation for
				      * obvious reasons, so you have
				      * to overload this function.
				      */
    virtual void solve_primal_problem () = 0;

				     /**
				      * This function is called by the
				      * manager object when solving the
				      * dual problem on this time level
				      * is needed. It is called after
				      * the #wake_up# function was
				      * called and before the #sleep#
				      * function will be called. There
				      * is a default implementation
				      * doing plain nothing since some
				      * problems may not need solving a
				      * dual problem. However, it
				      * will abort the program when
				      * being called anyway, since then you
				      * should really overload the function.
				      */
    virtual void solve_dual_problem ();

				     /**
				      * This function is called by the
				      * manager object when postprocessing
				      * this time level
				      * is needed. It is called after
				      * the #wake_up# function was
				      * called and before the #sleep#
				      * function will be called. There
				      * is a default implementation
				      * doing plain nothing since some
				      * problems may not need doing a
				      * postprocess step, e.g. if everything
				      * was already done when solving the
				      * primal problem. However, it
				      * will abort the program when
				      * being called anyway, since then you
				      * should really overload the function.
				      */
    virtual void postprocess_timestep ();
    
				     /**
				      * Compute the time difference to the
				      * last time step. If this timestep is
				      * the first one, this function will
				      * result in an exception. Though this
				      * behaviour seems a bit drastic, it
				      * is appropriate in most cases since
				      * if there is no previous time step
				      * you will need special treatment
				      * anyway.
				      */
    double get_backward_timestep () const;

				     /**
				      * Return the time difference to the next
				      * time step. With regard to the case
				      * that there is no next time step,
				      * the same applies as for the function
				      * above.
				      */
    double get_forward_timestep () const;
    
				     /**
				      * Exception
				      */
    DeclException0 (ExcGridNotDeleted);

				     /**
				      * Exception
				      */
    DeclException0 (ExcInternalError);
				     /**
				      * Exception
				      */
    DeclException0 (ExcPureVirtualFunctionCalled);
				     /**
				      * Exception
				      */
    DeclException0 (ExcCantComputeTimestep);
    
  protected:
				     /**				      
				      * Pointer to the previous time step object
				      * in the list.
				      */
    const TimeStepBase *previous_timestep;

				     /**				      
				      * Pointer to the next time step object
				      * in the list.
				      */
    const TimeStepBase *next_timestep;

				     /**
				      * Number of the sweep we are presently
				      * in. This number is reset by the time
				      * level manager before a sweep is
				      * started.
				      */
    unsigned int sweep_no;

				     /**
				      * Number of the time step, counted from
				      * zero onwards. This number is reset at
				      * the start of each sweep by the time
				      * level manager, since some time steps
				      * may have been inserted or deleted
				      * after the previous sweep.
				      */
    unsigned int timestep_no;

				     /**
				      * Discrete time this level operates on.
				      */
    const double time;

				     /**
				      * Variable storing whether
				      *  the solution * of a primal or
				      *  a dual problem is * actual,
				      *  or any of the other actions
				      *  specified in the enum. This
				      *  variable is set by the
				      *  #init_for_*# functions.
				      */
    NextAction next_action;
    
  private:
				     /**
				      * Reset the pointer to the previous time
				      * step; shall only be called by the
				      * time level manager object.
				      *
				      * This function is called at the set-up
				      * of the manager object and whenever
				      * a timestep is inserted or deleted.
				      */
    void set_previous_timestep (const TimeStepBase *previous);

				     /**
				      * Reset the pointer to the next time
				      * step; shall only be called by the
				      * time level manager object.
				      *
				      * This function is called at the set-up
				      * of the manager object and whenever
				      * a timestep is inserted or deleted.
				      */
    void set_next_timestep (const TimeStepBase *next);

				     /**
				      * Set the number this time step
				      * has in the list of timesteps.
				      * This function is called by the
				      * time step management object at
				      * the beginning of each sweep, to
				      * update information which may have
				      * changed due to addition or deleltion
				      * of time levels.
				      */
    void set_timestep_no (const unsigned int step_no);

				     /**
				      * Set the number of the sweep we are
				      * presently in. This function is
				      * called by the time level management
				      * object at start-up time of each
				      * sweep.
				      */
    void set_sweep_no (const unsigned int sweep_no);
    


				     // make the manager object a friend
    friend class TimeDependent;
};

	




/**
 * Specialization of #TimeStepBase# which addresses some aspects of grid handling.
 * In particular, this class is though to make handling of grids available that
 * are adaptively refined on each time step separately or with a loose coupling
 * between time steps. It also takes care of deleting and rebuilding grids when
 * memory resources are a point, through the #sleep# and #wake_up# functions
 * declared in the base class.
 *
 * @author Wolfgang Bangerth, 1999
 */
template <int dim>
class TimeStepBase_Tria :  public TimeStepBase
{
  public:
				     // forward declaration
    struct Flags;

				     /**
				      * Default constructor. Does nothing but
				      * throws an exception. We need to have
				      * such a constructor in order to satisfy
				      * the needs of derived classes, which
				      * take this class as a virtual base class
				      * and don't need to call this constructor
				      * of they are not terminal classes. The
				      * compiler would like to know a
				      * constructor to call anyway since it
				      * can't know that the class is not
				      * terminal.
				      */
    TimeStepBase_Tria ();
    
				     /**
				      * Constructor. Takes a coarse
				      * grid from which the grids on this
				      * time level will be derived and
				      * some flags steering the behaviour
				      * of this object.
				      *
				      * The ownership of the coarse grid
				      * stays with the creator of this
				      * object. However, it is locked
				      * from destruction to guarantee
				      * the lifetime of the coarse grid
				      * is longer than it is needed by
				      * this object.
				      */
    TimeStepBase_Tria (const double              time,
		       const Triangulation<dim> &coarse_grid,
		       const Flags              &flags);

				     /**
				      * Destructor. At present, this does
				      * not more than releasing the lock on
				      * the coarse grid triangulation given
				      * to the constructor.
				      */
    virtual ~TimeStepBase_Tria ();

    				     /**
				      * Reconstruct all the data that is
				      * needed for this time level to work.
				      * This function serves to reget all
				      * the variables and data structures
				      * to work again after they have been
				      * send to sleep some time before, or
				      * at the first time we visit this time
				      * level. In particular, it is used
				      * to reconstruct the triangulation,
				      * degree of freedom handlers, to reload
				      * data vectors in case they have been
				      * stored to disk, etc. By default,
				      * this function rebuilds the triangulation
				      * if the respective flag has been set to
				      * destroy it in the #sleep# function.
				      * It does so also the first time we
				      * hit this function and #wakeup_level#
				      * equals #flags.wakeup_level_to_build_grid#,
				      * independently of the value of the
				      * mentioned flag. (Actually, it does so
				      * whenever the triangulation pointer
				      * equals the Null pointer and the
				      * value of #wakeup_level# is right.)
				      *
				      * Since this is an important task, you
				      * should call this function from your
				      * own function, should you choose to
				      * overload it in your own class (which
				      * likely is the case), preferably at
				      * the beginning so that your function
				      * can take effect of the triangulation
				      * already existing.
				      */
    virtual void wake_up (const unsigned wakeup_level);

				     /**
				      * This is the opposite function
				      * to #wake_up#. It is used to
				      * delete data or save it to disk
				      * after they are no more needed
				      * for the present sweep. Typical
				      * kinds of data for this are
				      * data vectors, degree of
				      * freedom handlers,
				      * triangulation objects,
				      * etc. which occupy large
				      * amounts of memory and may
				      * therefore be externalized.
				      *
				      * By default, if the user specified so
				      * in the flags for this object, the
				      * triangulation is deleted and the
				      * refinement history saved such that
				      * the respective #wake_up# function can
				      * rebuild it. You should therefore call
				      * this function from your overloaded
				      * version, preferrably at the end so
				      * that your function can use the
				      * triangulation as long as ou need it.
				      */
    virtual void sleep (const unsigned);

  protected:
    
				     /**
				      * Triangulation used at this
				      * time level. Since this is
				      * something that every time
				      * stepping scheme needs to have,
				      * we can safely put it into the
				      * base class. Note that the
				      * triangulation is frequently
				      * deleted and rebuilt by the
				      * functions #sleep# and
				      * #wake_up# to save memory, if
				      * such a behaviour is specified
				      * in the #flags# structure.
				      */
    Triangulation<dim>       *tria;

				     /**
				      * Pointer to a grid which is to
				      * be used as the coarse grid for
				      * this time level.  This pointer
				      * is set through the
				      * constructor; ownership remains
				      * with the owner of this
				      * management object.  */
    const Triangulation<dim> &coarse_grid;
    
				     /**
				      * Some flags about how this time level
				      * shall behave. See the documentation
				      * of this struct to find out more about
				      * that.
				      */
    const Flags flags;

  private:
				     /**
				      * Vectors holding the refinement and
				      * coarsening flags of the different
				      * sweeps on this time level. The
				      * vectors therefore hold the history
				      * of the grid.
				      */
    vector<vector<bool> >   refine_flags, coarsen_flags;

				     /**
				      * The refinement
				      * flags of the triangulation are stored
				      * in a local variable thus allowing
				      * a restoration. The coarsening flags
				      * are also stored.
				      */
    void save_refine_flags ();

    				     /**
				      * Restore the grid according to the saved
				      * data. For this, the coarse grid is
				      * copied and the grid is stepwise
				      * rebuilt using the saved flags.
				      */
    void restore_grid ();
};



template <int dim>
struct TimeStepBase_Tria<dim>::Flags
{
				     /**
				      * Default constructor; yields an exception,
				      * so is not really usable.
				      */
    Flags ();
    
				     /**
				      * Constructor; see the different
				      * fields for a description of the
				      * meaning of the parameters.
				      */
    Flags (const unsigned int max_refinement_level,
	   const bool         delete_and_rebuild_tria,
	   const unsigned int wakeup_level_to_build_grid,
	   const unsigned int sleep_level_to_delete_grid);
    
				     /**
				      * Maximum level of a cell in the
				      * triangulation of a time level. If it
				      * is set to zero, then no limit is imposed
				      * on the number of refinements a coarse
				      * grid cell may undergo. Usually, this
				      * field is used, if for some reason you
				      * want to limit refinement in an
				      * adaptive process, for example to avoid
				      * overly large numbers of cells or to
				      * compare with grids which have a certain
				      * number of refinements.
				      */
    const unsigned int max_refinement_level;

				     /**
				      * This flag determines whether the
				      * #sleep# and #wake_up# functions shall
				      * delete and rebuild the triangulation.
				      * While for small problems, this is
				      * not necessary, for large problems
				      * it is indispensable to save memory.
				      * The reason for this is that there
				      * may be several hundred time levels
				      * in memory, each with its own
				      * triangulation, which may require
				      * large amounts if there are many
				      * cells on each. Having a total
				      * of 100.000.000 cells on all time
				      * levels taken together is not
				      * uncommon, which makes this flag
				      * understandable.
				      */
    const bool delete_and_rebuild_tria;

				     /**
				      * This number denotes the parameter to
				      * the #wake_up# function at which it
				      * shall rebuild the grid. Obviously,
				      * it shall be less than or equal to the
				      * #look_ahead# number passed to the
				      * time step management object; if it
				      * is equal, then the grid is rebuilt
				      * the first time the #wake_up# function
				      * is called. If #delete_and_rebuild_tria#
				      * is #false#, this number has no meaning.
				      */
    const unsigned int wakeup_level_to_build_grid;

				     /**
				      * This is the opposite flag to the one
				      * above: it determines at which call to
				      * #sleep# the grid shall be deleted.
				      */
    const unsigned int sleep_level_to_delete_grid;

				     /**
				      * Exception
				      */
    DeclException1 (ExcInvalidParameter,
		    int,
		    << "The parameter " << arg1 << " has an invalid value.");
				     /**
				      * Exception
				      */
    DeclException0 (ExcInternalError);
};




/*----------------------------   time-dependent.h     ---------------------------*/
/* end of #ifndef __time_dependent_H */
#endif
/*----------------------------   time-dependent.h     ---------------------------*/
