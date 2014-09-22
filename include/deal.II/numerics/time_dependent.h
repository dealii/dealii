// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#ifndef __deal2__time_dependent_h
#define __deal2__time_dependent_h


/*----------------------------   time-dependent.h     ---------------------------*/


#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/base/smartpointer.h>

#include <vector>
#include <utility>

DEAL_II_NAMESPACE_OPEN

// forward declarations
class TimeStepBase;
template <typename number> class Vector;
template <int dim, int spacedim> class Triangulation;

/**
 * This class provides an abstract interface to time dependent problems in that
 * it addresses some of the most annoying aspects of this class of problems:
 * data management. These problems frequently need large amounts of computer
 * ressources, most notably computing time, main memory and disk space.
 * Main memory reduction is often the most pressing need, methods to implement
 * it are almost always quite messy, though, quickly leading to code that
 * stores and reloads data at places scattered all over the program, and
 * which becomes unmaintainable sometimes. The present class tries to offer
 * a more structured interface, albeit simple, which emerged in my mind after
 * messing with my wave equation simulation for several months.
 *
 * The design of this class is mostly tailored for the solution of time
 * dependent partial differential equations where the computational
 * meshes may differ between each two timesteps and where the computations
 * on each time step take a rather long time compared with the overhead
 * of this class. Since no reference to the class of problems is made within
 * this class, it is not restricted to PDEs, though, and it seems likely that
 * a solver for large ordinary matrix differential equations may successfully
 * use the same setup and therefore this class.
 *
 *
 * <h3>Overview</h3>
 *
 * The general structure of a time dependent problem solver using a timestepping
 * scheme is about the following: we have a collection of time step objects
 * on which we solve our problem subsequently. In order to do so, we need
 * knowledge of the data on zero or several previous timesteps (when using single
 * or multiple step methods, that is) and maybe also some data of time steps
 * ahead (for example the computational grid on these). Dependening on the
 * problem in question, a second loop over all timesteps may be done solving
 * a dual problem, where the loop may run forward (one dual problem for each
 * time step) or backward (using a global dual problem). Within one of these
 * loops or using a separate loop, error estimators may be computed and the
 * grids may be refined. Each of these loops are initiated by a call preparing
 * each timestep object for the next loop, before actually starting the loop
 * itself.
 *
 * We will denote a complete set of all these loops with the term "sweep".
 * Since this library is mostly about adaptive methods, it is likely that the
 * last loop within a sweep will generate refined meshes and that we will
 * perform another sweep on these refined meshes. A total run will therefore
 * often be a sequence of several sweeps. The global setup therefore looks
 * like this:
 * @verbatim
 *    for sweep=0 to n_sweeps-1
 *    {
 *      for i=0 to n_timesteps-1
 *        initialize timestep i for this sweep, e.g. for setting up
 *        data structures, creating temporary files, etc.
 *
 *      for i=0 to n_timesteps-1
 *        prepare timestep i for loop 0
 *      for i=0 to n_timesteps-1
 *        perform loop 0 on timestep i   (e.g. solve primal problem)
 *
 *      for i=0 to n_timesteps-1
 *        prepare timestep i for loop 1
 *      for i=0 to n_timesteps-1
 *        perform loop 1 on timestep i   (e.g. solve dual problem)
 *
 *      for i=0 to n_timesteps-1
 *        prepare timestep i for loop 2
 *      for i=0 to n_timesteps-1
 *        perform loop 2 on timestep i   (e.g. compute error information)
 *
 *      ...
 *
 *      for i=0 to n_timesteps-1
 *        notify timestep i of the end of the sweep, e.g. for cleanups,
 *        deletion of temporary files, etc.
 *    }
 * @endverbatim
 * The user may specify that a loop shall run forward or backward (the latter
 * being needed for the solution of global dual problems, for example).
 *
 * Going from the global overview to a more local viewpoint, we note that when
 * a loop visits one timestep (e.g. to solve the primal or dual problem, or
 * to compute error information), we need information on this, one or more
 * previous time steps and zero or more timesteps in the future. However,
 * often it is not needed to know all information from these timesteps and
 * it is often a computational requirement to delete data at the first
 * possible time when it is no more needed. Likewise, data should be reloaded
 * at the latest time possible.
 *
 * In order to facilitate these principles, the concept of waking up and
 * letting sleep a time step object was developed. Assume we have a time
 * stepping scheme which needs to look ahead one time step and needs the
 * data of the last two time steps, the following pseudocode described
 * what the centeral loop function of this class will do when we move
 * from timestep @p n-1 to timestep @p n:
 * @verbatim
 *   wake up timestep n+1 with signal 1
 *   wake up timestep n with signal 0
 *   do computation on timestep n
 *   let timestep n sleep with signal 0
 *   let timestep n-1 sleep with signal 1
 *   let timestep n-2 sleep with signal 2
 *
 *   move from n to n+1
 * @endverbatim
 * The signal number here denotes the distance of the timestep being sent
 * the signal to the timestep where computations are done on. The calls to
 * the @p wake_up and @p sleep functions with signal 0 could in principle
 * be absorbed into the function doing the computation; we use these
 * redundant signals, however, in order to separate computations and data
 * management from each other, allowing to put all stuff around grid
 * management, data reload and storage into one set of functions and
 * computations into another.
 *
 * In the example above, possible actions might be: timestep <tt>n+1</tt> rebuilds
 * the computational grid (there is a specialized class which can do this
 * for you); timestep @p n builds matrices sets solution vectors to the right
 * size, maybe using an initial guess; then it does the computations; then
 * it deletes the matrices since they are not needed by subsequent timesteps;
 * timestep @p n-1 deletes those data vectors which are only needed by one
 * timestep ahead; timestep @p n-2 deletes the remaining vectors and deletes
 * the computational grid, somewhere storing information how to rebuild it
 * eventually.
 *
 * From the given sketch above, it is clear that each time step object sees
 * the following sequence of events:
 * @verbatim
 *   wake up with signal 1
 *   wake up signal 0
 *   do computation
 *   sleep with signal 0
 *   sleep with signal 1
 *   sleep with signal 2
 * @endverbatim
 * This pattern is repeated for each loop in each sweep.
 *
 * For the different loops within each sweep, the numbers of timesteps
 * to look ahead (i.e. the maximum signal number to the @p wake_up function)
 * and the look-behind (i.e. the maximum signal number to the @p sleep
 * function) can be chosen separately. For example, it is usually only
 * needed to look one time step behind when computing error estimation
 * (in some cases, it may vene be possible to not look ahead or back
 * at all, in which case only signals zero will be sent), while one
 * needs a look back of at least one for a timestepping method.
 *
 * Finally, a note on the direction of look-ahead and look-back is in
 * place: look-ahead always refers to the direction the loop is running
 * in, i.e. for loops running forward, @p wake_up is called for timestep
 * objects with a greater time value than the one previously computed on,
 * while @p sleep is called for timesteps with a lower time. If the loop
 * runs in the opposite direction, e.g. when solving a global dual
 * problem, this order is reversed.
 *
 *
 * <h3>Implementation</h3>
 *
 * The main loop of a program using this class will usually look like
 * the following one, taken modified from an application program that
 * isn't distributed as part of the library:
 * @code
 *   template <int dim>
 *   void TimeDependent_Wave<dim>::run_sweep (const unsigned int sweep_no)
 *   {
 *     start_sweep (sweep_no);
 *
 *     solve_primal_problem ();
 *
 *     if (compute_dual_problem)
 *       solve_dual_problem ();
 *
 *     postprocess ();
 *
 *     if (sweep_no != number_of_sweeps-1)
 *       refine_grids ();
 *
 *     write_statistics ();
 *
 *     end_sweep ();
 *   };
 *
 *
 *
 *   template <int dim>
 *   void WaveProblem<dim>::run ()
 *   {
 *     for (unsigned int sweep=0; sweep<number_of_sweeps; ++sweep)
 *       timestep_manager.run_sweep (sweep);
 *   };
 * @endcode
 * Here, @p timestep_manager is an object of type TimeDependent_Wave(), which
 * is a class derived from TimeDependent. @p start_sweep,
 * @p solve_primal_problem, @p solve_dual_problem, @p postprocess and @p end_sweep
 * are functions inherited from this class. They all do a loop over all
 * timesteps within this object and call the respective function on each of
 * these objects. For example, here are two of the functions as they are
 * implemented by the library:
 * @code
 *   void TimeDependent::start_sweep (const unsigned int s)
 *   {
 *     sweep_no = s;
 *
 *                                 // reset the number each
 *                                 // time step has, since some time
 *                                 // steps might have been added since
 *                                 // the last time we visited them
 *                                 //
 *                                 // also set the sweep we will
 *                                 // process in the sequel
 *     for (unsigned int step=0; step<timesteps.size(); ++step)
 *       {
 *         timesteps[step]->set_timestep_no (step);
 *         timesteps[step]->set_sweep_no (sweep_no);
 *       };
 *
 *     for (unsigned int step=0; step<timesteps.size(); ++step)
 *       timesteps[step]->start_sweep ();
 *   };
 *
 *
 *   void
 *   TimeDependent::solve_primal_problem ()
 *   {
 *     do_loop (mem_fun(&TimeStepBase::init_for_primal_problem),
 *              mem_fun(&TimeStepBase::solve_primal_problem),
 *              timestepping_data_primal,
 *              forward);
 *   };
 * @endcode
 * The latter function shows rather clear how most of the loops are
 * invoked (@p solve_primal_problem, @p solve_dual_problem, @p postprocess,
 * @p refine_grids and @p write_statistics all have this form, where the
 * latter two give functions of the derived timestep class, rather than
 * from the base class). The function TimeStepBase::init_for_primal_problem
 * and the respective ones for the other operations defined by that class
 * are only used to store the type of operation which the loop presently
 * performed will do.
 *
 * As can be seen, most of the work is done by the @p do_loop function of
 * this class, which takes the addresses of two functions which are used
 * to initialize all timestep objects for the loop and to actually perform
 * some action. The next parameter gives some information on the look-ahead
 * and look-back and the last one denotes in which direction the loop is
 * to be run.
 *
 * Using function pointers through the @p mem_fun functions provided by
 * the <tt>C++</tt> standard library, it is possible to do neat tricks, like
 * the following, also taken from the wave program, in this case from
 * the function @p refine_grids:
 * @code
 *   ...
 *   compute the thresholds for refinement
 *   ...
 *
 *   do_loop (mem_fun (&TimeStepBase_Tria<dim>::init_for_refinement),
 *            bind2nd (mem_fun1 (&TimeStepBase_Wave<dim>::refine_grid),
 *                     TimeStepBase_Tria<dim>::RefinementData (top_threshold,
 *                                                             bottom_threshold)),
 *            TimeDependent::TimeSteppingData (0,1),
 *            TimeDependent::forward);
 * @endcode
 * TimeStepBase_Wave::refine_grid is a function taking an argument, unlike
 * all the other functions used above within the loops. However, in this special
 * case the parameter was the same for all timesteps and known before the loop
 * was started, so we fixed it and made a function object which to the outside
 * world does not take parameters.
 *
 * Since it is the central function of this class, we finally present a
 * stripped down version of the @p do_loop method, which is shown in order
 * to provide a better understanding of the internals of this class. For
 * brevity we have omitted the parts that deal with backward running loops
 * as well as the checks whether wake-up and sleep operations act on timesteps
 * outside <tt>0..n_timesteps-1</tt>.
 * @code
 *   template <typename InitFunctionObject, typename LoopFunctionObject>
 *   void TimeDependent::do_loop (InitFunctionObject      init_function,
 *                           LoopFunctionObject      loop_function,
 *                           const TimeSteppingData &timestepping_data,
 *                           const Direction         direction)
 *   {
 *                                 // initialize the time steps for
 *                                 // a round of this loop
 *     for (unsigned int step=0; step<n_timesteps; ++step)
 *       init_function (static_cast<typename InitFunctionObject::argument_type>
 *                 (timesteps[step]));
 *
 *                                 // wake up the first few time levels
 *     for (int step=-timestepping_data.look_ahead; step<0; ++step)
 *       for (int look_ahead=0; look_ahead<=timestepping_data.look_ahead; ++look_ahead)
 *         timesteps[step+look_ahead]->wake_up(look_ahead);
 *
 *
 *     for (unsigned int step=0; step<n_timesteps; ++step)
 *       {
 *                                     // first thing: wake up the
 *                                     // timesteps ahead as necessary
 *         for (unsigned int look_ahead=0;
 *         look_ahead<=timestepping_data.look_ahead; ++look_ahead)
 *      timesteps[step+look_ahead]->wake_up(look_ahead);
 *
 *
 *                                     // actually do the work
 *         loop_function (static_cast<typename LoopFunctionObject::argument_type>
 *                   (timesteps[step]));
 *
 *                                     // let the timesteps behind sleep
 *         for (unsigned int look_back=0; look_back<=timestepping_data.look_back; ++look_back)
 *      timesteps[step-look_back]->sleep(look_back);
 *       };
 *
 *                                 // make the last few timesteps sleep
 *     for (int step=n_timesteps; n_timesteps+timestepping_data.look_back; ++step)
 *       for (int look_back=0; look_back<=timestepping_data.look_back; ++look_back)
 *         timesteps[step-look_back]->sleep(look_back);
 *   };
 * @endcode
 *
 *
 * @author Wolfgang Bangerth, 1999
 */
class TimeDependent
{
public:
  /**
   * Structure holding the two basic
   * entities that control a loop over
   * all time steps: how many time steps
   * ahead of the present one we shall
   * start waking up timestep objects
   * and how many timesteps behind
   * we shall call their @p sleep
   * method.
   */
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
     * the @p wake_up function for each
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
     * call the @p sleep function for
     * each time step.
     */
    const unsigned int look_back;
  };

  /**
   * Enum offering the different directions
   * in which a loop executed by
   * @p do_loop may be run.
   */
  enum Direction
  {
    forward, backward
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
   * given to the <tt>insert_*</tt> and
   * @p add_timestep functions, i.e.
   * it will delete the objects doing
   * the computations on each time step.
   */
  virtual ~TimeDependent ();

  /**
   * Add a timestep at any position. The
   * position is a pointer to an existing
   * time step object, or a null pointer
   * denoting the end of the timestep
   * sequence. If @p position is non-null,
   * the new time step will be inserted
   * before the respective element.
   *
   * Note that by giving an object
   * to this function, the
   * TimeDependent object assumes
   * ownership of the object; it will
   * therefore also take care of
   * deletion of the objects its manages.
   *
   * There is another function,
   * @p add_timestep, which inserts a
   * time step at the end of the list.
   *
   * Note that this function does not
   * change the timestep numbers stored
   * within the other timestep objects,
   * nor does it set the timestep number
   * of this new timestep. This is only
   * done upon calling the @p start_sweep
   * function. In not changing the timestep
   * numbers, it is simpler to operate
   * on a space-time triangulation since
   * one can always use the timestep numbers
   * that were used in the previous sweep.
   */
  void insert_timestep (const TimeStepBase *position,
                        TimeStepBase       *new_timestep);

  /**
   * Just like @p insert_timestep, but
   * insert at the end.
   *
   * This mechanism usually will result
   * in a set-up loop like this
   * @code
   * for (i=0; i<N; ++i)
   *   manager.add_timestep(new MyTimeStep());
   * @endcode
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
   *
   * Note that this function does
   * not change the timestep
   * numbers stored within the
   * other timestep objects. This
   * is only done upon calling the
   * @p start_sweep function. In not
   * changing the timestep numbers,
   * it is simpler to operate on a
   * space-time triangulation since
   * one can always use the
   * timestep numbers that were
   * used in the previous sweep.
   */
  void delete_timestep (const unsigned int position);

  /**
   * Solve the primal problem; uses the
   * functions @p init_for_primal_problem
   * and @p solve_primal_problem of the
   * TimeStepBase class through the
   * @p do_loop function of this class.
   *
   * Look ahead and look back are
   * determined by the @p timestepping_data_primal
   * object given to the constructor.
   */
  void solve_primal_problem ();

  /**
   * Solve the dual problem; uses the
   * functions @p init_for_dual_problem
   * and @p solve_dual_problem of the
   * TimeStepBase class through the
   * @p do_loop function of this class.
   *
   * Look ahead and look back are
   * determined by the @p timestepping_data_dual
   * object given to the constructor.
   */
  void solve_dual_problem ();

  /**
   * Do a postprocessing round; uses the
   * functions @p init_for_postprocessing
   * and @p postprocess of the
   * TimeStepBase class through the
   * @p do_loop function of this class.
   *
   * Look ahead and look back are
   * determined by the @p timestepping_data_postprocess
   * object given to the constructor.
   */
  void postprocess ();

  /**
   * Do a loop over all timesteps, call the
   * @p init_function at the beginning and
   * the @p loop_function of each time step.
   * The @p timestepping_data determine how
   * many timesteps in front and behind
   * the present one the @p wake_up and
   * @p sleep functions are called.
   *
   * To see how this function work, note that
   * the function @p solve_primal_problem only
   * consists of a call to
   * <tt>do_loop (mem_fun(&TimeStepBase::init_for_primal_problem),
   *    mem_fun(&TimeStepBase::solve_primal_problem),
   *    timestepping_data_primal, forward);</tt>.
   *
   * Note also, that the given class from which
   * the two functions are taken needs not
   * necessarily be TimeStepBase, but it
   * could also be a derived class, that is
   * @p static_castable from a TimeStepBase.
   * The function may be a virtual function
   * (even a pure one) of that class, which
   * should help if the actual class where it
   * is implemented is one which is derived
   * through virtual base classes and thus
   * unreachable by @p static_cast from the
   * TimeStepBase class.
   *
   * Instead of using the above form, you can
   * equally well use
   * <tt>bind2nd(mem_fun1(&X::unary_function), arg)</tt>
   * which lets the @p do_loop
   * function call the given function with
   * the specified parameter. Note that you
   * need to bind the second parameter since
   * the first one implicitly contains
   * the object which the function is to
   * be called for.
   */
  template <typename InitFunctionObject, typename LoopFunctionObject>
  void do_loop (InitFunctionObject      init_function,
                LoopFunctionObject      loop_function,
                const TimeSteppingData &timestepping_data,
                const Direction         direction);


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
   * It also calls the @p start_sweep
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
   *
   * The default implementation of this
   * function calls @p start_sweep on all
   * time step objects.
   */
  virtual void start_sweep (const unsigned int sweep_no);

  /**
   * Analogon to the above
   * function, calling @p end_sweep
   * of each time step object. The
   * same applies with respect to
   * the @p virtualness of this
   * function as for the previous
   * one.
   *
   * @note This function does not
   * guarantee that @p end_sweep is
   * called for successive time
   * steps successively, rather the order of
   * time step objects for which the
   * function is called is
   * arbitrary. You should
   * therefore not assume that that
   * function has been called for
   * previous time steps
   * already. If in multithread
   * mode, the @p end_sweep function
   * of several time steps may be
   * called at once, so you should
   * use synchronization
   * mechanisms if your program
   * requires so.
   */
  virtual void end_sweep ();

  /**
   * @deprecated Use the function without an argument.
   * @arg n_threads This argument is ignored.
   */
  virtual void end_sweep (const unsigned int n_threads) DEAL_II_DEPRECATED;

  /**
   * Determine an estimate for the
   * memory consumption (in bytes)
   * of this object.
   */
  std::size_t memory_consumption () const;

  /**
   * Exception.
   */
  DeclException0 (ExcInvalidPosition);

protected:
  /**
   * Vector holding pointers to the time
   * level objects. This is the main data
   * this object operates on. Note that
   * this object takes possession of the
   * objects pointed to by the pointers
   * in this collection.
   */
  std::vector<SmartPointer<TimeStepBase,TimeDependent> > timesteps;

  /**
   * Number of the present sweep. This is
   * reset by the @p start_sweep function
   * called at the outset of each sweep.
   */
  unsigned int sweep_no;

  /**
   * Some flags telling the
   * @p solve_primal_problem function what to
   * do. See the documentation of this struct
   * for more information.
   */
  const TimeSteppingData timestepping_data_primal;

  /**
   * Some flags telling the
   * @p solve_dual_problem function what to
   * do. See the documentation of this struct
   * for more information.
   */
  const TimeSteppingData timestepping_data_dual;

  /**
   * Some flags telling the
   * @p postprocess function what to
   * do. See the documentation of this struct
   * for more information.
   */
  const TimeSteppingData timestepping_data_postprocess;

private:

  /**
   * Do the work of <tt>end_sweep()</tt>
   * for some timesteps only. This
   * is useful in multithread mode.
   */
  void end_sweep (const unsigned int begin_timestep,
                  const unsigned int end_timestep);
};



/**
 * Base class for a time step in time dependent problems. This class provides
 * barely more than the basic framework, defining the necessary virtual
 * functions (namely @p sleep and @p wake_up), the interface to previous
 * and following grids, and some functions to be called before a new loop
 * over all time steps is started.
 *
 * @author Wolfgang Bangerth, 1999
 */
class TimeStepBase : public Subscriptor
{
public:
  /**
   * Enum denoting the type of problem
   * which will have to be solved next.
   */
  enum SolutionState
  {
    primal_problem = 0x0,
    dual_problem   = 0x1,
    postprocess    = 0x2
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
  virtual void wake_up (const unsigned int);

  /**
   * This is the opposite function
   * to @p wake_up. It is used to
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
  virtual void sleep (const unsigned int);

  /**
   * This function is called each time
   * before a new sweep is started. You
   * may want to set up some fields needed
   * in the course of the computations,
   * and so on. You should take good care,
   * however, not to install large objects,
   * which should be deferred until the
   * @p wake_up function is called.
   *
   * A typical action of this function
   * would be sorting out names of
   * temporary files needed in the
   * process of solving, etc.
   *
   * At the time this function is called,
   * the values of @p timestep_no, @p sweep_no
   * and the pointer to the previous and
   * next time step object already have
   * their correct value.
   *
   * The default implementation of this
   * function does nothing.
   */
  virtual void start_sweep ();

  /**
   * This is the analogon to the above
   * function, but it is called at the
   * end of a sweep. You will usually want
   * to do clean-ups in this function,
   * such as deleting temporary files
   * and the like.
   */
  virtual void end_sweep ();

  /**
   * Before the primal problem is
   * solved on each time level, this
   * function is called (i.e. before the
   * solution takes place on the first
   * time level). By default, this function
   * sets the @p next_action variable of
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
   * the @p wake_up function was
   * called and before the @p sleep
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
   * the @p wake_up function was
   * called and before the @p sleep
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
   * the @p wake_up function was
   * called and before the @p sleep
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
   * Return the time value of this time
   * step.
   */
  double get_time () const;

  /**
   * Return the number of this time
   * step. Note that this number may vary
   * between different sweeps, if timesteps
   * are added or deleted.
   */
  unsigned int get_timestep_no () const;

  /**
   * Compute the time difference to the
   * last time step. If this timestep is
   * the first one, this function will
   * result in an exception. Though this
   * behaviour seems a bit drastic, it
   * is appropriate in most cases since
   * if there is no previous time step
   * you will need special treatment
   * anyway and this way no invalid
   * value is returned which could lead
   * to wrong but unnoticed results of
   * your computation. (The only sensible
   * value to return in that case would
   * not be zero, since valid computation
   * can be done with that, but would
   * be a denormalized value such as @p NaN.
   * However, there is not much difference
   * in finding that the results of a
   * computation are all denormalized values
   * or in getting an exception; in the
   * latter case you at least get the exact
   * place where your problem lies.)
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
   * Determine an estimate for the
   * memory consumption (in bytes)
   * of this object.
   *
   * You will want to overload
   * this function in derived
   * classes to compute the
   * amount memory used by the
   * derived class, and add the
   * result of this function to
   * your result.
   */
  virtual std::size_t memory_consumption () const;

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
   * Variable storing whether the
   * solution of a primal or a dual
   * problem is actual, or any of
   * the other actions
   * specified. This variable is
   * set by the <tt>init_for_*</tt>
   * functions.
   */
  unsigned int next_action;

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


  /**
   * Copy constructor. I can see no reason
   * why someone might want to use it, so
   * I don't provide it. Since this class
   * has pointer members, making it private
   * prevents the compiler to provide it's
   * own, incorrect one if anyone chose to
   * copy such an object.
   */
  TimeStepBase (const TimeStepBase &);

  /**
   * Copy operator. I can see no reason
   * why someone might want to use it, so
   * I don't provide it. Since this class
   * has pointer members, making it private
   * prevents the compiler to provide it's
   * own, incorrect one if anyone chose to
   * copy such an object.
   */
  TimeStepBase &operator = (const TimeStepBase &);

  // make the manager object a friend
  friend class TimeDependent;
};




/**
 * Namespace in which some classes are declared that encapsulate flags
 * for the TimeStepBase_Tria() class. These used to be local data
 * types of that class, but some compilers choked on some aspects, so
 * we put them into a namespace of their own.
 *
 * @author Wolfgang Bangerth, 2001
 */
namespace TimeStepBase_Tria_Flags
{
  /**
   * This structure is used to tell the TimeStepBase_Tria() class how grids should
   * be handled. It has flags defining the moments where grids shall be
   * re-made and when they may be deleted. Also, one variable states whether
   * grids should be kept in memory or should be deleted between to uses to
   * save memory.
   */
  template <int dim>
  struct Flags
  {
    /**
     * Default constructor; yields
     * an exception, so is not
     * really usable.
     */
    Flags ();

    /**
     * Constructor; see the different
     * fields for a description of the
     * meaning of the parameters.
     */
    Flags (const bool         delete_and_rebuild_tria,
           const unsigned int wakeup_level_to_build_grid,
           const unsigned int sleep_level_to_delete_grid);

    /**
     * This flag determines whether
     * the @p sleep and
     * @p wake_up functions shall
     * delete and rebuild the
     * triangulation.  While for
     * small problems, this is not
     * necessary, for large
     * problems it is indispensable
     * to save memory.  The reason
     * for this is that there may
     * be several hundred time
     * levels in memory, each with
     * its own triangulation, which
     * may require large amounts if
     * there are many cells on
     * each. Having a total of
     * 100.000.000 cells on all
     * time levels taken together
     * is not uncommon, which makes
     * this flag understandable.
     */
    const bool delete_and_rebuild_tria;

    /**
     * This number denotes the
     * parameter to the @p wake_up
     * function at which it shall
     * rebuild the grid. Obviously,
     * it shall be less than or
     * equal to the @p look_ahead
     * number passed to the time
     * step management object; if
     * it is equal, then the grid
     * is rebuilt the first time
     * the @p wake_up function is
     * called. If
     * @p delete_and_rebuild_tria
     * is @p false, this number
     * has no meaning.
     */
    const unsigned int wakeup_level_to_build_grid;

    /**
     * This is the opposite flag to
     * the one above: it determines
     * at which call to * @p sleep
     * the grid shall be deleted.
     */
    const unsigned int sleep_level_to_delete_grid;

    /**
     * Exception
     */
    DeclException1 (ExcInvalidParameter,
                    int,
                    << "The parameter " << arg1 << " has an invalid value.");
  };



  /**
   * This structure is used to tell the TimeStepBase_Tria() class how grids should
   * be refined. Before we explain all the different variables, fist some terminology:
   * <ul>
   * <li> Correction: after having flagged some cells of the triangulation for
   *   following some given criterion, we may want to change the number of flagged
   *   cells on this grid according to another criterion that the number of cells
   *   may be only a certain fraction more or less then the number of cells on
   *   the previous grid. This change of refinement flags will be called
   *   "correction" in the sequel.
   * <li> Adaption: in order to make the change between one grid and the next not
   *   to large, we may want to flag some additional cells on one of the two
   *   grids such that there are not too grave differences. This process will
   *   be called "adaption".
   * </ul>
   *
   *
   * <h3>Description of flags</h3>
   *
   * <ul>
   * <li> @p max_refinement_level: Cut the refinement of cells at a given level.
   *   This flag does not influence the flagging of cells, so not more cells
   *   on the coarser levels are flagged than usual. Rather, the flags are all
   *   set, but when it comes to the actual refinement, the maximum refinement
   *   level is truncated.
   *
   *   This option is only really useful when you want to compare global
   *   refinement with adaptive refinement when you don't want the latter
   *   to refine more than the global refinement.
   *
   * <li> @p first_sweep_with_correction: When using cell number correction
   *   as defined above, it may be worth while to start with this only in
   *   later sweeps, not already in the first one. If this variable is
   *   zero, then start with the first sweep, else with a higher one. The
   *   rationale for only starting later is that we do not want to block the
   *   development of grids at the beginning and only impose restrictions in
   *   the sweeps where we start to be interested in the actual results of
   *   the computations.
   *
   * <li> @p min_cells_for_correction: If we want a more free process of
   *   grid development, we may want to impose less rules for grids with few
   *   cells also. This variable sets a lower bound for the cell number of
   *   grids where corrections are to be performed.
   *
   * <li> @p cell_number_corridor_top: Fraction of the number of cells by
   *   which the number of cells of one grid may be higher than that on the
   *   previous grid. Common values are 10 per cent (i.e. 0.1). The naming
   *   of the variable results from the goal to define a target corridor
   *   for the number of cells after refinement has taken place.
   *
   * <li> @p cell_number_corridor_bottom: Fraction of the number of cells by
   *   which the number of cells of one grid may be lower than that on the
   *   previous grid. Common values are 5 per cent (i.e. 0.05). Usually this
   *   number will be smaller than @p cell_number_corridor_top since an
   *   increase of the number of cells is not harmful (though it increases
   *   the numerical amount of work needed to solve the problem) while a
   *   sharp decrease may reduce the accuracy of the final result even if
   *   the time steps computed before the decrease were computed to high
   *   accuracy.
   *
   *   Note however, that if you compute the dual problem as well, then the time
   *   direction is reversed, so the two values defining the cell number
   *   corridor should be about equal.
   *
   * <li> @p correction_relaxations: This is a list of pairs of number with the
   *   following meaning: just as for @p min_cells_for_correction, it may be
   *   worth while to reduce the requirements upon grids if the have few cells.
   *   The present variable stores a list of cell numbers along with some values
   *   which tell us that the cell number corridor should be enlarged by a
   *   certain factor. For example, if this list was <tt>((100 5) (200 3) (500 2))</tt>,
   *   this would mean that for grids with a cell number below 100, the
   *   <tt>cell_number_corridor_*</tt> variables are to be multiplied by 5 before they
   *   are applied, for cell numbers below 200 they are to be multiplied by 3,
   *   and so on.
   *
   *   @p correction_relaxations is actually a vector of such list. Each entry
   *   in this vector denotes the relaxation rules for one sweep. The last
   *   entry defines the relaxation rules for all following sweeps. This
   *   scheme is adopted to allow for stricter corrections in later sweeps
   *   while the relaxations may be more generous in the first sweeps.
   *
   *   There is a static variable @p default_correction_relaxations which you
   *   can use as a default value. It is an empty list and thus defines no
   *   relaxations.
   *
   * <li> @p cell_number_correction_steps: Usually, if you want the number of
   *   cells to be corrected, the target corridor for the cell number is computed
   *   and some additional cells are flagged or flags are removed. But since
   *   the cell number resulting after flagging and deflagging can not be
   *   easily computed, it will usually not be within the corridor. We therefore
   *   need to iteratively get to our goal. Usually, three or four iterations are
   *   needed, but using this variable, you can reduce the allowed number of
   *   iterations; breaking the loop after two iterations yields good results
   *   regularly. Setting the variable to zero will result in no correction
   *   steps at all.
   *
   * <li> @p mirror_flags_to_previous_grid: If a cell on the present grid is
   *   flagged for refinement, also flag the corresponding cell on the previous
   *   grid. This is useful if, for example, error indicators are computed for
   *   space-time cells, but are stored for the second grid only. Now, since the
   *   first grid has the same contributions to the indicators as the second, it
   *   may be useful to flag both if necessary. This is done if the present
   *   variable is set.
   *
   * <li> @p adapt_grids: adapt the present grid to the previous one in the sense
   *   defined above. What is actually done here is the following: if going from
   *   the previous to the present grid would result in double refinement or
   *   double coarsening of some cells, then we try to flag these cells for
   *   refinement or coarsening such as to avoid the double step. Obviously, more
   *   than double refinement of coarsening is also caught.
   *
   *   Grid adaption can try to avoid such changes between two grids, but it can
   *   never promise that they don't occur. This is because the next grid may
   *   change the present one, but then again there may be jumps in refinement
   *   level between the present and the previous one; this could only be avoided
   *   by looping iteratively through all grids, back and forth, until nothing
   *   changes anymore, which is obviously impossible if there are many time steps
   *   with very large grids.
   * </ul>
   */
  template <int dim>
  struct RefinementFlags
  {
    /**
     * Typedef of a data type
     * describing some relaxations
     * of the correction process.
     * See the general description
     * of this class for more
     * information.
     */
    typedef std::vector<std::vector<std::pair<unsigned int, double> > >   CorrectionRelaxations;

    /**
     * Default values for the relaxations:
     * no relaxations.
     */
    static CorrectionRelaxations default_correction_relaxations;

    /**
     * Constructor. The default
     * values are chosen such that
     * almost no restriction on the
     * mesh refinement is imposed.
     */
    RefinementFlags (const unsigned int max_refinement_level         = 0,
                     const unsigned int first_sweep_with_correction  = 0,
                     const unsigned int min_cells_for_correction     = 0,
                     const double cell_number_corridor_top           = (1<<dim),
                     const double cell_number_corridor_bottom        = 1,
                     const CorrectionRelaxations &correction_relaxations = CorrectionRelaxations(),
                     const unsigned int cell_number_correction_steps = 0,
                     const bool mirror_flags_to_previous_grid        = false,
                     const bool adapt_grids                          = false);

    /**
     * Maximum level of a cell in
     * the triangulation of a time
     * level. If it is set to zero,
     * then no limit is imposed on
     * the number of refinements a
     * coarse grid cell may
     * undergo. Usually, this field
     * is used, if for some reason
     * you want to limit refinement
     * in an adaptive process, for
     * example to avoid overly
     * large numbers of cells or to
     * compare with grids which
     * have a certain number of
     * refinements.
     */
    const unsigned int  max_refinement_level;

    /**
     * First sweep to perform cell
     * number correction steps on;
     * for sweeps before, cells are
     * only flagged and no
     * number-correction to
     * previous grids is performed.
     */
    const unsigned int  first_sweep_with_correction;


    /**
     * Apply cell number correction
     * with the previous time level
     * only if there are more than
     * this number of cells.
     */
    const unsigned int  min_cells_for_correction;

    /**
     * Fraction by which the number
     * of cells on a time level may
     * differ from the number on
     * the previous time level
     * (first: top deviation,
     * second: bottom deviation).
     */
    const double        cell_number_corridor_top;

    /**
     * @ref cell_number_corridor_top
     */
    const double        cell_number_corridor_bottom;

    /**
     * List of relaxations to the
     * correction step.
     */
    const std::vector<std::vector<std::pair<unsigned int,double> > > correction_relaxations;

    /**
     * Number of iterations to be
     * performed to adjust the
     * number of cells on a time
     * level to those on the
     * previous one. Zero means: do
     * no such iteration.
     */
    const unsigned int  cell_number_correction_steps;

    /**
     * Flag all cells which are
     * flagged on this timestep for
     * refinement on the previous
     * one also. This is useful in
     * case the error indicator was
     * computed by integration over
     * time-space cells, but are
     * now associated to a grid on
     * a discrete time level. Since
     * the error contribution comes
     * from both grids, however, it
     * is appropriate to refine
     * both grids.
     *
     * Since the previous grid does
     * not mirror the flags to the
     * one before it, this does not
     * lead to an almost infinite
     * growth of cell numbers. You
     * should use this flag with
     * cell number correction
     * switched on only, however.
     *
     * Mirroring is done after cell
     * number correction is done,
     * but before grid adaption, so
     * the cell number on this grid
     * is not noticeably influenced
     * by the cells flagged
     * additionally on the previous
     * grid.
     */
    const bool          mirror_flags_to_previous_grid;

    /**
     * Adapt this grid to the
     * previous one.
     */
    const bool          adapt_grids;

    /**
     * Exception
     */
    DeclException1 (ExcInvalidValue,
                    double,
                    << "The following value does not fulfill the requirements: " << arg1);
  };



  /**
   * Structure given to the actual refinement function, telling it which
   * thresholds to take for coarsening and refinement. The actual refinement
   * criteria are loaded by calling the virtual function
   * @p get_tria_refinement_criteria.
   */
  template <int dim>
  struct RefinementData
  {
    /**
     * Constructor
     */
    RefinementData (const double         refinement_threshold,
                    const double         coarsening_threshold=0);

    /**
     * Threshold for refinement:
     * cells having a larger value
     * will be refined (at least in
     * the first round; subsequent
     * steps of the refinement
     * process may flag other cells
     * as well or remove the flag
     * from cells with a criterion
     * higher than this threshold).
     */
    const double         refinement_threshold;

    /**
     * Same threshold for
     * coarsening: cells with a
     * smaller threshold will be
     * coarsened if possible.
     */
    const double         coarsening_threshold;

    /**
     * Exception
     */
    DeclException1 (ExcInvalidValue,
                    double,
                    << "The following value does not fulfill the requirements: " << arg1);
  };
}




/**
 * Specialisation of TimeStepBase which addresses some aspects of grid handling.
 * In particular, this class is thought to make handling of grids available that
 * are adaptively refined on each time step separately or with a loose coupling
 * between time steps. It also takes care of deleting and rebuilding grids when
 * memory resources are a point, through the @p sleep and @p wake_up functions
 * declared in the base class.
 *
 * In addition to that, it offers functions which do some rather hairy refinement
 * rules for time dependent problems, trying to avoid too much change in the grids
 * between subsequent time levels, while also trying to retain the freedom of
 * refining each grid separately. There are lots of flags and numbers controlling
 * this function, which might drastically change the behaviour of the function -- see
 * the description of the flags for further information.
 *
 * @author Wolfgang Bangerth, 1999; large parts taken from the wave program, by Wolfgang Bangerth 1998
 */
template <int dim>
class TimeStepBase_Tria : public TimeStepBase
{
public:
  /**
   * Typedef the data types of the
   * TimeStepBase_Tria_Flags()
   * namespace into local scope.
   */
  typedef typename TimeStepBase_Tria_Flags::Flags<dim>           Flags;
  typedef typename TimeStepBase_Tria_Flags::RefinementFlags<dim> RefinementFlags;
  typedef typename TimeStepBase_Tria_Flags::RefinementData<dim>  RefinementData;


  /**
   * Extension of the enum in the base
   * class denoting the next action to be
   * done.
   */
  enum SolutionState
  {
    grid_refinement = 0x1000
  };


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
   *
   * You need to give the general flags
   * structure to this function since it
   * is needed anyway; the refinement
   * flags can be omitted if you do
   * not intend to call the refinement
   * function of this class.
   */
  TimeStepBase_Tria (const double              time,
                     const Triangulation<dim, dim> &coarse_grid,
                     const Flags              &flags,
                     const RefinementFlags    &refinement_flags = RefinementFlags());

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
   * destroy it in the @p sleep function.
   * It does so also the first time we
   * hit this function and @p wakeup_level
   * equals <tt>flags.wakeup_level_to_build_grid</tt>,
   * independently of the value of the
   * mentioned flag. (Actually, it does so
   * whenever the triangulation pointer
   * equals the Null pointer and the
   * value of @p wakeup_level is right.)
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
  virtual void wake_up (const unsigned int wakeup_level);

  /**
   * This is the opposite function
   * to @p wake_up. It is used to
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
   * the respective @p wake_up function can
   * rebuild it. You should therefore call
   * this function from your overloaded
   * version, preferably at the end so
   * that your function can use the
   * triangulation as long as ou need it.
   */
  virtual void sleep (const unsigned int);

  /**
   * Do the refinement according to the
   * flags passed to the constructor of this
   * object and the data passed to this
   * function. For a description of the
   * working of this function refer to the
   * general documentation of this class.
   *
   * In fact, this function does not
   * actually refine or coarsen the
   * triangulation, but only sets the
   * respective flags. This is done because
   * usually you will not need the grid
   * immediately afterwards but only
   * in the next sweep, so it suffices
   * to store the flags and rebuild it
   * the next time we need it. Also, it
   * may be that the next time step
   * would like to add or delete some
   * flags, so we have to wait anyway
   * with the use of this grid.
   */
  void refine_grid (const RefinementData data);

  /**
   * Respective init function for the
   * refinement loop; does nothing in
   * the default implementation, apart from
   * setting @p next_action to
   * @p grid_refinement but can
   * be overloaded.
   */
  virtual void init_for_refinement ();

  /**
   * Virtual function that should fill
   * the vector with the refinement
   * criteria for the present triangulation.
   * It is used within the @p refine_grid
   * function to get the criteria for
   * the present time step, since they
   * can't be passed through its
   * argument when using the loop of
   * the time step management object.
   */
  virtual void get_tria_refinement_criteria (Vector<float> &criteria) const = 0;

  /**
   * The refinement
   * flags of the triangulation are stored
   * in a local variable thus allowing
   * a restoration. The coarsening flags
   * are also stored.
   */
  void save_refine_flags ();

  /**
   * Determine an estimate for the
   * memory consumption (in bytes)
   * of this object.
   *
   * You will want to overload
   * this function in derived
   * classes to compute the
   * amount memory used by the
   * derived class, and add the
   * result of this function to
   * your result.
   */
  virtual std::size_t memory_consumption () const;

  /**
   * Exception
   */
  DeclException0 (ExcGridNotDeleted);

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
   * functions @p sleep and
   * @p wake_up to save memory, if
   * such a behaviour is specified
   * in the @p flags structure.
   */
  SmartPointer<Triangulation<dim, dim>,TimeStepBase_Tria<dim> > tria;

  /**
   * Pointer to a grid which is to
   * be used as the coarse grid for
   * this time level.  This pointer
   * is set through the
   * constructor; ownership remains
   * with the owner of this
   * management object.
   */
  SmartPointer<const Triangulation<dim, dim>,TimeStepBase_Tria<dim> > coarse_grid;

  /**
   * Some flags about how this time level
   * shall behave. See the documentation
   * of this struct to find out more about
   * that.
   */
  const Flags flags;

  /**
   * Flags controlling the refinement
   * process; see the documentation of the
   * respective structure for more
   * information.
   */
  const RefinementFlags refinement_flags;

private:
  /**
   * Vectors holding the refinement and
   * coarsening flags of the different
   * sweeps on this time level. The
   * vectors therefore hold the history
   * of the grid.
   */
  std::vector<std::vector<bool> >   refine_flags;

  /**
   * @ref refine_flags
   */
  std::vector<std::vector<bool> >   coarsen_flags;

  /**
   * Restore the grid according to the saved
   * data. For this, the coarse grid is
   * copied and the grid is stepwise
   * rebuilt using the saved flags.
   */
  void restore_grid ();
};





/*----------------------------- template functions ------------------------------*/

template <typename InitFunctionObject, typename LoopFunctionObject>
void TimeDependent::do_loop (InitFunctionObject      init_function,
                             LoopFunctionObject      loop_function,
                             const TimeSteppingData &timestepping_data,
                             const Direction         direction)
{
  // the following functions looks quite
  // disrupted due to the recurring switches
  // for forward and backward running loops.
  //
  // I chose to switch at every place where
  // it is needed, since it is so easy
  // to overlook something when you change
  // some code at one place when it needs
  // to be changed at a second place, here
  // for the other direction, also.

  const unsigned int n_timesteps = timesteps.size();

  // initialize the time steps for
  // a round of this loop
  for (unsigned int step=0; step<n_timesteps; ++step)
    switch (direction)
      {
      case forward:
        init_function (static_cast<typename InitFunctionObject::argument_type>
                       (&*timesteps[step]));
        break;
      case backward:
        init_function (static_cast<typename InitFunctionObject::argument_type>
                       (&*timesteps[n_timesteps-step-1]));
        break;
      };


  // wake up the first few time levels
  for (int step=-timestepping_data.look_ahead; step<0; ++step)
    for (int look_ahead=0;
         look_ahead<=static_cast<int>(timestepping_data.look_ahead); ++look_ahead)
      switch (direction)
        {
        case forward:
          if (step+look_ahead >= 0)
            timesteps[step+look_ahead]->wake_up(look_ahead);
          break;
        case backward:
          if (n_timesteps-(step+look_ahead) < n_timesteps)
            timesteps[n_timesteps-(step+look_ahead)]->wake_up(look_ahead);
          break;
        };


  for (unsigned int step=0; step<n_timesteps; ++step)
    {
      // first thing: wake up the
      // timesteps ahead as necessary
      for (unsigned int look_ahead=0;
           look_ahead<=timestepping_data.look_ahead; ++look_ahead)
        switch (direction)
          {
          case forward:
            if (step+look_ahead < n_timesteps)
              timesteps[step+look_ahead]->wake_up(look_ahead);
            break;
          case backward:
            if (n_timesteps > (step+look_ahead))
              timesteps[n_timesteps-(step+look_ahead)-1]->wake_up(look_ahead);
            break;
          };


      // actually do the work
      switch (direction)
        {
        case forward:
          loop_function (static_cast<typename LoopFunctionObject::argument_type>
                         (&*timesteps[step]));
          break;
        case backward:
          loop_function (static_cast<typename LoopFunctionObject::argument_type>
                         (&*timesteps[n_timesteps-step-1]));
          break;
        };

      // let the timesteps behind sleep
      for (unsigned int look_back=0;
           look_back<=timestepping_data.look_back; ++look_back)
        switch (direction)
          {
          case forward:
            if (step>=look_back)
              timesteps[step-look_back]->sleep(look_back);
            break;
          case backward:
            if (n_timesteps-(step-look_back) <= n_timesteps)
              timesteps[n_timesteps-(step-look_back)-1]->sleep(look_back);
            break;
          };
    };

  // make the last few timesteps sleep
  for (int step=n_timesteps;
       step<static_cast<int>(n_timesteps+timestepping_data.look_back); ++step)
    for (int look_back=0;
         look_back<=static_cast<int>(timestepping_data.look_back); ++look_back)
      switch (direction)
        {
        case forward:
          if ((step-look_back >= 0)
              &&
              (step-look_back < static_cast<int>(n_timesteps)))
            timesteps[step-look_back]->sleep(look_back);
          break;
        case backward:
          if ((step-look_back >= 0)
              &&
              (step-look_back < static_cast<int>(n_timesteps)))
            timesteps[n_timesteps-(step-look_back)-1]->sleep(look_back);
          break;
        };
}

DEAL_II_NAMESPACE_CLOSE

/*----------------------------   time-dependent.h     ---------------------------*/
#endif
/*----------------------------   time-dependent.h     ---------------------------*/
