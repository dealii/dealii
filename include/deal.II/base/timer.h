// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2017 by the deal.II authors
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

#ifndef dealii_timer_h
#define dealii_timer_h

#include <deal.II/base/config.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/thread_management.h>
#include <deal.II/base/utilities.h>

#include <chrono>
#include <list>
#include <map>
#include <string>

DEAL_II_NAMESPACE_OPEN

/**
 * A clock, compatible with the <code>std::chrono</code> notion of a clock,
 * whose now() method returns a time point indicating the amount of CPU time
 * that the current process has used.
 */
struct CPUClock
{
  /**
   * Duration type. Windows measures CPU times, by default, in multiples of
   * 1/64th of a second and POSIX uses microseconds, so go with microseconds
   * for uniformity.
   */
  typedef std::chrono::microseconds duration;

  /**
   * Signed integral type used to store the value returned by count().
   */
  typedef duration::rep rep;

  /**
   * Ratio representing the length of a period (in seconds).
   */
  typedef duration::period period;

  /**
   * Time point type.
   */
  typedef std::chrono::time_point<CPUClock, duration> time_point;

  /**
   * Boolean indicating that the clock monotonically increases.
   */
  static const bool is_steady = true;

  /**
   * Return the amount of CPU time that the current process has
   * used. Unfortunately, this requires platform-specific calls, so this
   * function returns 0 on platforms that are neither Windows nor POSIX.
   */
  static time_point now() noexcept;
};

/**
 * The Timer class provides a way to measure both the amount of wall time
 * (i.e., the amount of time elapsed on a wall clock) and the amount of CPU
 * time that certain sections of an application have used. This class also
 * offers facilities for synchronizing the elapsed time across an MPI
 * communicator.
 *
 * <h3>Usage</h3>
 *
 * The Timer class can be started and stopped several times. It stores both
 * the amount of time elapsed over the last start-stop cycle, or <em>lap</em>,
 * as well as the total time elapsed over all laps. Here is an example:
 *
 * @code
 *   Timer timer; // creating a timer also starts it
 *
 *   // do some complicated computations here
 *   // ...
 *
 *   timer.stop();
 *
 *   std::cout << "Elapsed CPU time: " << timer.cpu_time() << " seconds.\n";
 *   std::cout << "Elapsed wall time: " << timer.wall_time() << " seconds.\n";
 *
 *   // reset timer for the next thing it shall do
 *   timer.reset();
 * @endcode
 *
 * Alternatively, you can also restart the timer instead of resetting it. The
 * times between successive calls to start() and stop() (i.e., the laps) will
 * then be accumulated. The usage of this class is also explained in the
 * step-28 tutorial program.
 *
 * @note The TimerOutput (combined with TimerOutput::Scope) class provide a
 * convenient way to time multiple named sections and summarize the output.
 *
 * @note Implementation of this class is system dependent. In particular, CPU
 * times are accumulated from summing across all threads and will usually
 * exceed the wall times.
 *
 * @ingroup utilities
 * @author G. Kanschat, W. Bangerth, M. Kronbichler, D. Wells
 */
class Timer
{
public:
  /**
   * Constructor. Sets the accumulated times at zero and calls Timer::start().
   */
  Timer ();

  /**
   * Constructor specifying that CPU times should be summed over the given
   * communicator. If @p sync_lap_times is <code>true</code> then the Timer
   * will set the elapsed wall and CPU times over the last lap to their
   * maximum values across the provided communicator. This synchronization is
   * only performed if Timer::stop() is called before the timer is queried for
   * time duration values.
   *
   * This constructor calls Timer::start().
   *
   * @note The timer is stopped before the synchronization over the
   * communicator occurs; the extra cost of the synchronization is not
   * measured.
   */
  Timer (MPI_Comm mpi_communicator,
         const bool sync_lap_times = false);

  /**
   * Return a reference to the data structure with global timing information
   * for the last lap. This structure does not contain meaningful values until
   * Timer::stop() has been called.
   *
   * @deprecated Use Timer::get_last_lap_wall_time_data() instead, which
   * returns a reference to the same structure.
   */
  DEAL_II_DEPRECATED
  const Utilities::MPI::MinMaxAvg &get_data() const;

  /**
   * Return a reference to the data structure containing basic statistics on
   * the last lap's wall time measured across all MPI processes in the given
   * communicator. This structure does not contain meaningful values until
   * Timer::stop() has been called.
   */
  const Utilities::MPI::MinMaxAvg &get_last_lap_wall_time_data() const;

  /**
   * Return a reference to the data structure containing basic statistics on
   * the accumulated wall time measured across all MPI processes in the given
   * communicator. This structure does not contain meaningful values until
   * Timer::stop() has been called.
   *
   * @deprecated Use Timer::get_accumulated_wall_time_data() instead, which
   * returns a reference the same structure.
   */
  DEAL_II_DEPRECATED
  const Utilities::MPI::MinMaxAvg &get_total_data() const;

  /**
   * Return a reference to the data structure containing basic statistics on
   * the accumulated wall time measured across all MPI processes in the given
   * communicator. This structure does not contain meaningful values until
   * Timer::stop() has been called.
   */
  const Utilities::MPI::MinMaxAvg &get_accumulated_wall_time_data() const;

  /**
   * Prints the data returned by get_data(), i.e. for the last lap,
   * to the given stream.
   *
   * @deprecated Use Timer::print_last_lap_wall_time_data() instead, which
   * prints the same information.
   */
  template <class StreamType>
  DEAL_II_DEPRECATED
  void print_data(StreamType &stream) const;

  /**
   * Print the data returned by Timer::get_last_lap_wall_time_data() to the
   * given stream.
   */
  template <class StreamType>
  void print_last_lap_wall_time_data(StreamType &stream) const;

  /**
   * Prints the data returned by get_total_data(), i.e. for the total run,
   * to the given stream.
   *
   * @deprecated Use Timer::print_accumulated_wall_time_data() instead, which
   * prints the same information.
   */
  template <class StreamType>
  DEAL_II_DEPRECATED
  void print_total_data(StreamType &stream) const;

  /**
   * Print the data returned by Timer::get_accumulated_wall_time_data() to the
   * given stream.
   */
  template <class StreamType>
  void print_accumulated_wall_time_data(StreamType &stream) const;

  /**
   * Begin measuring a new lap. If <code>sync_lap_times</code> is
   * <code>true</code> then an MPI barrier is used to ensure that all
   * processes begin the lap at the same wall time.
   */
  void start ();

  /**
   * Stop the timer. This updates the lap times and accumulated times. If
   * <code>sync_lap_times</code> is <code>true</code> then the lap times are
   * synchronized over all processors in the communicator (i.e., the lap times
   * are set to the maximum lap time).
   *
   * Return the accumulated CPU time in seconds.
   */
  double stop ();

  /**
   * Stop the timer, if it is running, and reset all measured values to their
   * default states.
   */
  void reset ();

  /**
   * Equivalent to calling Timer::reset() followed by calling Timer::start().
   */
  void restart();

  /**
   * Access to the current CPU time without stopping the timer. The elapsed
   * time is returned in units of seconds.
   *
   * @deprecated Use cpu_time() instead.
   */
  DEAL_II_DEPRECATED
  double operator() () const;

  /**
   * Return the current accumulated wall time (including the current lap, if
   * the timer is running) in seconds without stopping the timer.
   */
  double wall_time () const;

  /**
   * Return the wall time of the last lap in seconds. The timer is not stopped
   * by this function.
   */
  double last_wall_time() const;

  /**
   * Return the accumulated CPU time (including the current lap, if the timer
   * is running) in seconds without stopping the timer.
   *
   * If an MPI communicator is provided to the constructor then the returned
   * value is the sum of all accumulated CPU times over all processors in the
   * communicator.
   */
  double cpu_time() const;

  /**
   * Return the CPU time of the last lap in seconds. The timer is not stopped
   * by this function.
   */
  double last_cpu_time() const;

  /**
   * Return the wall time taken between the last start()/stop() call.
   *
   * @deprecated Use last_wall_time() instead.
   */
  DEAL_II_DEPRECATED
  double get_lap_time () const;

private:
  /**
   * The Timer class stores timing information for two different clocks: a
   * wall clock and a CPU usage clock. Since the logic for handling both
   * clocks is, in most places, identical, we collect the relevant
   * measurements for each clock into this <code>struct</code>.
   *
   * @tparam clock_type_ The type of the clock whose measurements are being
   * stored. This class should conform to the usual clock interface expected
   * by <code>std::chrono</code> (i.e., the correct <code>typedef</code>s and
   * a static <code>now()</code> method).
   */
  template <class clock_type_>
  struct ClockMeasurements
  {
    /**
     * Store the clock type.
     */
    typedef clock_type_ clock_type;

    /**
     * The time point type of the provided clock.
     */
    typedef typename clock_type::time_point time_point_type;

    /**
     * The duration type of the provided clock.
     */
    typedef typename clock_type::duration duration_type;

    /**
     * The time point corresponding to the start of the current lap. This is
     * obtained by calling <code>clock_type::now()</code>.
     */
    time_point_type current_lap_start_time;

    /**
     * The accumulated time over several laps.
     */
    duration_type accumulated_time;

    /**
     * The duration of the last lap.
     */
    duration_type last_lap_time;

    /**
     * Constructor. Sets <code>current_lap_start_time</code> to the current
     * clock time and the durations to zero.
     */
    ClockMeasurements();

    /**
     * Reset the clock by setting <code>current_lap_start_time</code> to the
     * current clock time and the durations to zero.
     */
    void reset();
  };

  /**
   * typedef for the wall clock.
   */
  typedef std::chrono::steady_clock wall_clock_type;

  /**
   * typedef for the CPU clock.
   */
  typedef CPUClock cpu_clock_type;

  /**
   * Collection of wall time measurements.
   */
  ClockMeasurements<wall_clock_type> wall_times;

  /**
   * Collection of CPU time measurements.
   */
  ClockMeasurements<cpu_clock_type> cpu_times;

  /**
   * Whether or not the timer is presently running.
   */
  bool running;

  /**
   * The communicator over which various time values are synchronized and
   * combined: see the documentation of the relevant constructor for
   * additional information.
   */
  MPI_Comm mpi_communicator;

  /**
   * Store whether or not the wall time and CPU time are synchronized across
   * the communicator in Timer::start() and Timer::stop().
   */
  bool sync_lap_times;

  /**
   * A structure for parallel wall time measurement that includes the minimum,
   * maximum, and average over all processors known to the MPI communicator of
   * the last lap time.
   */
  Utilities::MPI::MinMaxAvg last_lap_wall_time_data;

  /**
   * A structure for parallel wall time measurement that includes the minimum
   * time recorded among all processes, the maximum time as well as the
   * average time defined as the sum of all individual times divided by the
   * number of MPI processes in the MPI_Comm for the total run time.
   */
  Utilities::MPI::MinMaxAvg accumulated_wall_time_data;
};



//TODO: The following class is not thread-safe
/**
 * This class can be used to generate formatted output from time measurements
 * of different subsections in a program. It is possible to create several
 * sections that perform certain aspects of the program. A section can be
 * entered several times. By changing the options in OutputFrequency and
 * OutputType, the user can choose whether output should be generated every
 * time a section is joined or just in the end of the program. Moreover, it is
 * possible to show CPU times, wall times or both.
 *
 * <h3>Usage</h3>
 *
 * Use of this class could be as follows:
 * @code
 *   TimerOutput timer (std::cout, TimerOutput::summary,
 *                      TimerOutput::wall_times);
 *
 *   timer.enter_subsection ("Setup dof system");
 *   setup_dofs();
 *   timer.leave_subsection();
 *
 *   timer.enter_subsection ("Assemble");
 *   assemble_system_1();
 *   timer.leave_subsection();
 *
 *   timer.enter_subsection ("Solve");
 *   solve_system_1();
 *   timer.leave_subsection();
 *
 *   timer.enter_subsection ("Assemble");
 *   assemble_system_2();
 *   timer.leave_subsection();
 *
 *   timer.enter_subsection ("Solve");
 *   solve_system_2();
 *   timer.leave_subsection();
 *
 *   // do something else...
 * @endcode
 * When run, this program will return an output like this:
 * @code
 * +---------------------------------------------+------------+------------+
 * | Total wallclock time elapsed since start    |      88.8s |            |
 * |                                             |            |            |
 * | Section                         | no. calls |  wall time | % of total |
 * +---------------------------------+-----------+------------+------------+
 * | Assemble                        |         2 |      19.7s |        22% |
 * | Solve                           |         2 |      3.03s |       3.4% |
 * | Setup dof system                |         1 |      3.97s |       4.5% |
 * +---------------------------------+-----------+------------+------------+
 * @endcode
 * The output will see that we entered the assembly and solve section twice,
 * and reports how much time we spent there. Moreover, the class measures the
 * total time spent from start to termination of the TimerOutput object. In
 * this case, we did a lot of other stuff, so that the time proportions of the
 * functions we measured are far away from 100 percent.
 *
 *
 * <h3>Using scoped timers</h3>
 *
 * The scheme above where you have to have calls to
 * TimerOutput::enter_subsection() and TimerOutput::leave_subsection() is
 * awkward if the sections in between these calls contain <code>return</code>
 * statements or may throw exceptions. In that case, it is easy to forget that
 * one nevertheless needs to leave the section somehow, somewhere. An easier
 * approach is to use "scoped" sections. This is a variable that when you
 * create it enters a section, and leaves the section when you destroy it. If
 * this is a variable local to a particular scope (a code block between curly
 * braces) and you leave this scope due to a <code>return</code> statements or
 * an exception, then the variable is destroyed and the timed section is left
 * automatically. Consequently, we could have written the code piece above as
 * follows, with exactly the same result but now exception-safe:
 * @code
 *   TimerOutput timer (std::cout, TimerOutput::summary,
 *                      TimerOutput::wall_times);
 *
 *   {
 *     TimerOutput::Scope timer_section(timer, "Setup dof system");
 *     setup_dofs();
 *   }
 *
 *   {
 *     TimerOutput::Scope timer_section(timer, "Assemble");
 *     assemble_system_1();
 *   }
 *
 *   {
 *     TimerOutput::Scope timer_section(timer, "Solve");
 *     solve_system_1();
 *   }
 *
 *   {
 *     TimerOutput::Scope timer_section(timer, "Assemble");
 *     assemble_system_2();
 *   }
 *
 *   {
 *     TimerOutput::Scope timer_section(timer, "Solve");
 *     solve_system_2();
 *   }
 *
 *   // do something else...
 * @endcode
 *
 *
 * <h3>Usage in parallel programs using MPI</h3>
 *
 * In a parallel program built on MPI, using the class in a way such as the
 * one shown above would result in a situation where each process times the
 * corresponding sections and then outputs the resulting timing information at
 * the end. This is annoying since you'd get a lot of output -- one set of
 * timing information from each processor.
 *
 * This can be avoided by only letting one processor generate screen output,
 * typically by using an object of type ConditionalOStream instead of
 * <code>std::cout</code> to write to screen (see, for example, step-17,
 * step-18, step-32 and step-40, all of which use this method).
 *
 * This way, only a single processor outputs timing information, typically the
 * first process in the MPI universe. However, if you take the above code
 * snippet as an example, imagine what would happen if
 * <code>setup_dofs()</code> is fast on processor zero and slow on at least
 * one of the other processors; and if the first thing
 * <code>assemble_system_1()</code> does is something that requires all
 * processors to communicate. In this case, on processor zero, the timing
 * section with name <code>"Setup dof system"</code> will yield a short run
 * time on processor zero, whereas the section <code> "Assemble"</code> will
 * take a long time: not because <code>assemble_system_1()</code> takes a
 * particularly long time, but because on the processor on which we time (or,
 * rather, the one on which we generate output) happens to have to wait for a
 * long time till the other processor is finally done with
 * <code>setup_dofs()</code> and starts to participate in
 * <code>assemble_system_1()</code>. In other words, the timing that is
 * reported is unreliable because it reflects run times from other processors.
 * Furthermore, the run time of this section on processor zero has nothing to
 * do with the run time of the section on other processors but instead with
 * the run time of <i>the previous section</i> on another processor.
 *
 * The usual way to avoid this is to introduce a barrier into the parallel
 * code just before we start and stop timing sections. This ensures that all
 * processes are at the same place and the timing information then reflects
 * the maximal run time across all processors. To achieve this, you need to
 * initialize the TimerOutput object with an MPI communicator object, for
 * example as in the following code:
 * @code
 *   TimerOutput timer (MPI_COMM_WORLD,
 *                      pcout,
 *                      TimerOutput::summary,
 *                      TimerOutput::wall_times);
 * @endcode
 * Here, <code>pcout</code> is an object of type ConditionalOStream that makes
 * sure that we only generate output on a single processor. See the step-32,
 * step-40, and step-42 tutorial programs for this kind of usage of this class.
 *
 * @ingroup utilities
 * @author M. Kronbichler, 2009.
 */
class TimerOutput
{
public:
  /**
   * Helper class to enter/exit sections in TimerOutput be constructing a
   * simple scope-based object. The purpose of this class is explained in the
   * documentation of TimerOutput.
   */
  class Scope
  {
  public:
    /**
     * Enter the given section in the timer. Exit automatically when calling
     * stop() or destructor runs.
     */
    Scope(dealii::TimerOutput &timer_, const std::string &section_name);

    /**
     * Destructor calls stop()
     */
    ~Scope();

    /**
     * In case you want to exit the scope before the destructor is executed,
     * call this function.
     */
    void stop();

  private:
    /**
     * Reference to the TimerOutput object
     */
    dealii::TimerOutput &timer;

    /**
     * Name of the section we need to exit
     */
    const std::string section_name;

    /**
     * Do we still need to exit the section we are in?
     */
    bool in;
  };

  /**
   * An enumeration data type that describes whether to generate output every
   * time we exit a section, just in the end, both, or never.
   */
  enum OutputFrequency
  {
    /**
     * Generate output after every call.
     */
    every_call,
    /**
     * Generate output in summary at the end.
     */
    summary,
    /**
     * Generate output both after every call and in summary at the end.
     */
    every_call_and_summary,
    /**
     * Never generate any output.
     */
    never
  };

  /**
   * An enumeration data type that describes the type of data to return
   * when fetching the data from the timer.
   */
  enum OutputData
  {
    /**
     * Output CPU times.
     */
    total_cpu_time,
    /**
     * Output wall clock times.
     */
    total_wall_time,
    /**
     * Output number of calls.
     */
    n_calls
  };

  /**
   * An enumeration data type that describes whether to show CPU times, wall
   * times, or both CPU and wall times whenever we generate output.
   */
  enum OutputType
  {
    /**
     * Output CPU times.
     */
    cpu_times,
    /**
     * Output wall clock times.
     */
    wall_times,
    /**
     * Output both CPU and wall clock times.
     */
    cpu_and_wall_times
  };

  /**
   * Constructor.
   *
   * @param stream The stream (of type std::ostream) to which output is
   * written.
   * @param output_frequency A variable indicating when output is to be
   * written to the given stream.
   * @param output_type A variable indicating what kind of timing the output
   * should represent (CPU or wall time).
   */
  TimerOutput (std::ostream          &stream,
               const OutputFrequency  output_frequency,
               const OutputType       output_type);

  /**
   * Constructor.
   *
   * @param stream The stream (of type ConditionalOstream) to which output is
   * written.
   * @param output_frequency A variable indicating when output is to be
   * written to the given stream.
   * @param output_type A variable indicating what kind of timing the output
   * should represent (CPU or wall time).
   */
  TimerOutput (ConditionalOStream    &stream,
               const OutputFrequency  output_frequency,
               const OutputType       output_type);

  /**
   * Constructor that takes an MPI communicator as input. A timer constructed
   * this way will sum up the CPU times over all processors in the MPI network
   * for calculating the CPU time, or take the maximum over all processors,
   * depending on the value of @p output_type . See the documentation of this
   * class for the rationale for this constructor and an example.
   *
   * @param mpi_comm An MPI communicator across which we should accumulate or
   * otherwise synchronize the timing information we produce on every MPI
   * process.
   * @param stream The stream (of type std::ostream) to which output is
   * written.
   * @param output_frequency A variable indicating when output is to be
   * written to the given stream.
   * @param output_type A variable indicating what kind of timing the output
   * should represent (CPU or wall time). In this parallel context, when this
   * argument selects CPU time, then times are accumulated over all processes
   * participating in the MPI communicator. If this argument selects wall
   * time, then reported times are the maximum over all processors' run times
   * for this section. (The latter is computed by placing an
   * <code>MPI_Barrier</code> call before starting and stopping the timer for
   * each section.
   */
  TimerOutput (MPI_Comm               mpi_comm,
               std::ostream          &stream,
               const OutputFrequency  output_frequency,
               const OutputType       output_type);

  /**
   * Constructor that takes an MPI communicator as input. A timer constructed
   * this way will sum up the CPU times over all processors in the MPI network
   * for calculating the CPU time, or take the maximum over all processors,
   * depending on the value of @p output_type . See the documentation of this
   * class for the rationale for this constructor and an example.
   *
   * @param mpi_comm An MPI communicator across which we should accumulate or
   * otherwise synchronize the timing information we produce on every MPI
   * process.
   * @param stream The stream (of type ConditionalOstream) to which output is
   * written.
   * @param output_frequency A variable indicating when output is to be
   * written to the given stream.
   * @param output_type A variable indicating what kind of timing the output
   * should represent (CPU or wall time). In this parallel context, when this
   * argument selects CPU time, then times are accumulated over all processes
   * participating in the MPI communicator. If this argument selects wall
   * time, then reported times are the maximum over all processors' run times
   * for this section. (The latter is computed by placing an
   * <code>MPI_Barrier</code> call before starting and stopping the timer for
   * each section.)
   */
  TimerOutput (MPI_Comm               mpi_comm,
               ConditionalOStream    &stream,
               const OutputFrequency  output_frequency,
               const OutputType       output_type);

  /**
   * Destructor. Calls print_summary() in case the option for writing the
   * summary output is set.
   */
  ~TimerOutput();

  /**
   * Open a section by given a string name of it. In case the name already
   * exists, that section is entered once again and times are accumulated.
   */
  void enter_subsection (const std::string &section_name);

  /**
   * Same as @p enter_subsection.
   */
  void enter_section (const std::string &section_name);

  //TODO: make some of these functions DEPRECATED (I would keep enter/exit_section)

  /**
   * Leave a section. If no name is given, the last section that was entered
   * is left.
   */
  void leave_subsection (const std::string &section_name = std::string());

  /**
   * Same as @p leave_subsection.
   */
  void exit_section (const std::string &section_name = std::string());

  /**
   * Get a map with the collected data of the specified type for each subsection
   */
  std::map<std::string, double> get_summary_data (const OutputData kind) const;

  /**
   * Print a formatted table that summarizes the time consumed in the various
   * sections.
   */
  void print_summary () const;

  /**
   * By calling this function, all output can be disabled. This function
   * together with enable_output() can be useful if one wants to control the
   * output in a flexible way without putting a lot of <tt>if</tt> clauses in
   * the program.
   */
  void disable_output ();

  /**
   * This function re-enables output of this class if it was previously
   * disabled with disable_output(). This function together with
   * disable_output() can be useful if one wants to control the output in a
   * flexible way without putting a lot of <tt>if</tt> clauses in the program.
   */
  void enable_output ();

  /**
   * Resets the recorded timing information.
   */
  void reset ();

private:
  /**
   * When to output information to the output stream.
   */
  OutputFrequency output_frequency;

  /**
   * Whether to show CPU times, wall times, or both CPU and wall times.
   */
  OutputType output_type;


  /**
   * A timer object for the overall run time. If we are using MPI, this timer
   * also accumulates over all MPI processes.
   */
  Timer              timer_all;

  /**
   * A structure that groups all information that we collect about each of the
   * sections.
   */
  struct Section
  {
    Timer  timer;
    double total_cpu_time;
    double total_wall_time;
    unsigned int n_calls;
  };

  /**
   * A list of all the sections and their information.
   */
  std::map<std::string, Section> sections;

  /**
   * The stream object to which we are to output.
   */
  ConditionalOStream out_stream;

  /**
   * A boolean variable that sets whether output of this class is currently on
   * or off.
   */
  bool output_is_enabled;

  /**
   * A list of the sections that have been entered and not exited. The list is
   * kept in the order in which sections have been entered, but elements may
   * be removed in the middle if an argument is given to the exit_section()
   * function.
   */
  std::list<std::string> active_sections;

  /**
   * mpi communicator
   */
  MPI_Comm            mpi_communicator;

  /**
   * A lock that makes sure that this class gives reasonable results even when
   * used with several threads.
   */
  Threads::Mutex mutex;
};



/* ---------------- inline functions ----------------- */


inline
void Timer::restart ()
{
  reset();
  start();
}



inline
const Utilities::MPI::MinMaxAvg &
Timer::get_data() const
{
  return last_lap_wall_time_data;
}



inline
const Utilities::MPI::MinMaxAvg &
Timer::get_last_lap_wall_time_data() const
{
  return last_lap_wall_time_data;
}



inline
const Utilities::MPI::MinMaxAvg &
Timer::get_total_data() const
{
  return accumulated_wall_time_data;
}



inline
const Utilities::MPI::MinMaxAvg &
Timer::get_accumulated_wall_time_data() const
{
  return accumulated_wall_time_data;
}



template <class StreamType>
inline
void
Timer::print_data(StreamType &stream) const
{
  print_last_lap_wall_time_data(stream);
}



template <class StreamType>
inline
void
Timer::print_last_lap_wall_time_data(StreamType &stream) const
{
  const Utilities::MPI::MinMaxAvg &statistic = get_last_lap_wall_time_data();
  stream << statistic.max << " wall,"
         << " max @" << statistic.max_index << ", min=" << statistic.min << " @"
         << statistic.min_index << ", avg=" << statistic.avg << std::endl;
}



template <class StreamType>
inline
void
Timer::print_total_data(StreamType &stream) const
{
  print_accumulated_wall_time_data(stream);
}



template <class StreamType>
inline
void
Timer::print_accumulated_wall_time_data(StreamType &stream) const
{
  const Utilities::MPI::MinMaxAvg statistic = get_accumulated_wall_time_data();
  stream << statistic.max << " wall,"
         << " max @" << statistic.max_index << ", min=" << statistic.min << " @"
         << statistic.min_index << ", avg=" << statistic.avg << std::endl;
}



inline
void
TimerOutput::enter_section (const std::string &section_name)
{
  enter_subsection(section_name);
}



inline
void
TimerOutput::exit_section (const std::string &section_name)
{
  leave_subsection(section_name);
}

inline
TimerOutput::Scope::Scope(dealii::TimerOutput &timer_,
                          const std::string &section_name_)
  :
  timer(timer_),
  section_name(section_name_),
  in(true)
{
  timer.enter_section(section_name);
}



inline
void
TimerOutput::Scope::stop()
{
  if (!in) return;
  in=false;

  timer.exit_section(section_name);
}


DEAL_II_NAMESPACE_CLOSE

#endif
