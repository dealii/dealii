/*----------------------------   time-dependent.h     ---------------------------*/
/*      $Id$                 */
#ifndef __time_dependent_H
#define __time_dependent_H
/*----------------------------   time-dependent.h     ---------------------------*/


#include <base/exceptions.h>
#include <base/subscriptor.h>
#include <basic/forward-declarations.h>

#include <vector>



template <int dim>
class TimeDependent 
{
  public:
    
    struct TimeSteppingData;
    
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
    void insert_timestep (TimeStepBase<dim> *new_timestep,
			  const unsigned int position);

				     /**
				      * Just like #insert_timestep#, but
				      * insert at the end.
				      */
    void add_timestep (TimeStepBase<dim> *new_timestep);
    

    void solve_primal_problem ();

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
		    << "Can't insert time step at position " << arg1
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
    vector<TimeStepBase<dim>*> timesteps;

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
};




template <int dim>
struct TimeDependent<dim>::TimeSteppingData
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
				      */
    const unsigned int look_back;
};
    




template <int dim>
class TimeStepBase : public Subscriptor
{
  public:
				     // forward declaration
    struct Flags;

				     /**
				      * Enum denoting the type of problem
				      * which will have to be solved next.
				      */
    enum ProblemType {
	  primal_problem,
	  dual_problem
    };
    
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
    TimeStepBase (const Triangulation<dim> &coarse_grid,
		  const Flags              &flags);

				     /**
				      * Destructor. At present, this does
				      * not more than releasing the lock on
				      * the coarse grid triangulation given
				      * to the constructor.
				      */
    virtual ~TimeStepBase ();
    
    virtual void wake_up ();
    virtual void sleep ();

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
				      * sets the #problem_type# variable of
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
				      * Variable storing whether the solution
				      * of a primal or a dual problem is
				      * actual. This variable is set by the
				      * two functions
				      * #init_for_{primal,dual}_problem#.
				      */
    ProblemType problem_type;
    
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


				     // make the manager object a friend
    template <> friend class TimeDependent<dim>;
};

	

template <int dim>
struct TimeStepBase<dim>::Flags
{
				     /**
				      * Constructor; see the different
				      * fields for a description of the
				      * meaning of the parameters.
				      */
    Flags (const unsigned int max_refinement_level,
	   const bool delete_and_rebuild_tria);
    
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
};


/*----------------------------   time-dependent.h     ---------------------------*/
/* end of #ifndef __time_dependent_H */
#endif
/*----------------------------   time-dependent.h     ---------------------------*/
