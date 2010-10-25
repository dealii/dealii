//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2009, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__timer_h
#define __deal2__timer_h

#include <base/config.h>
#include <base/conditional_ostream.h>
#include <base/thread_management.h>

#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
#include <mpi.h>
#include <base/utilities.h>
#endif

#include <string>
#include <list>
#include <map>

DEAL_II_NAMESPACE_OPEN

/**
 * This is a very simple class which provides information about both the CPU
 * time and the wallclock time elapsed since the timer was started last
 * time. Information is retrieved from the system on the basis of clock cycles
 * since last time the computer was booted for the CPU time. The wall time is
 * based on the system clock accessed by @p gettimeofday, with a typical
 * accuracy of 0.01 ms on linux systems.
 *
 *
 * <h3>Usage</h3>
 *
 * Use of this class is as you might expect by looking at the member
 * functions:
 * @code
 *   Timer timer;
 *   timer.start ();
 *
 *   // do some complicated computations here
 *   ...
 *
 *   timer.stop ();
 *
 *   std::cout << "Elapsed CPU time: " << timer() << " seconds.";
 *   std::cout << "Elapsed wall time: " << timer.wall_time() << " seconds.";
 *
 *   // reset timer for the next thing it shall do
 *   timer.reset();
 * @endcode
 *
 * Alternatively, you can also restart the timer instead of resetting
 * it. The times between successive calls to start() / stop() will then be
 * accumulated. The usage of this class is also explained in the
 * step-12 and step-29 tutorial programs.
 *
 * @note Implementation of this class is system dependent. In case
 * multithreaded routines (matrix-vector products, error estimators, etc.) are
 * used, the CPU time is accumulated from all the children.
 *
 * @ingroup utilities
 * @author G. Kanschat, W. Bangerth, M. Kronbichler
 */
class Timer
{
  public:
				     /**
				      * Constructor. Starts the timer at 0 sec.
				      */
    Timer ();

#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
				     /**
				      * Constructor that takes an MPI
				      * communicator as input. A timer
				      * constructed this way will sum up the
				      * CPU times over all processors in the
				      * MPI network when requested by the
				      * operator ().
				      *
				      * Starts the timer at 0 sec.
				      *
				      * If @p sync_wall_time is true, the wall
				      * time is synchronized between all CPUs
				      * using a MPI_Barrier() and a collective
				      * operation. Note that this only works
				      * if you stop() the timer before
				      * querying for the wall time. The time
				      * for the MPI operations are not
				      * included in the timing but may slow
				      * down your program.
				      *
				      * This constructor is only available
				      * if the deal.II compiler is an MPI
				      * compiler.
				      */
    Timer (MPI_Comm mpi_communicator,
	   bool sync_wall_time = false);

				     /**
				      * Returns a reference to the data
				      * structure with global timing
				      * information. Filled after calling
				      * stop().
				      */
    const Utilities::System::MinMaxAvg & get_data() const;

				     /**
				      * Prints the data to the given stream.
				      */
    template <class STREAM>
    void print_data(STREAM & stream) const;


#endif

				     /**
				      * Re-start the timer at the point where
				      * it was stopped. This way a cumulative
				      * measurement of time is possible.
				      */
    void start ();

				     /**
				      * Sets the current time as next
				      * starting time and return the
				      * elapsed time in seconds.
				      */
    double stop ();

				     /**
				      * Stop the timer if necessary and reset
				      * the elapsed time to zero.
				      */
    void reset ();

				     /**
				      * Resets the elapsed time to zero and
				      * starts the timer. This corresponds to
				      * calling @p reset() and @p start() on
				      * the Timer object.
				      */
    void restart();

				     /**
				      * Access to the current CPU time
				      * without disturbing time
				      * measurement. The elapsed time is
				      * returned in units of seconds.
				      */
    double operator() () const;

				     /**
				      * Access to the current wall time
				      * without disturbing time
				      * measurement. The elapsed time is
				      * returned in units of seconds.
				      */
    double wall_time () const;

  private:

				     /**
				      * Value of the user time when start()
				      * was called the last time or when the
				      * object was created and no stop() was
				      * issued in between.
				      */
    double              start_time;


				     /**
				      * Similar to #start_time, but
				      * needed for children threads
				      * in multithread mode. Value of
				      * the user time when start()
				      * was called the last time or
				      * when the object was created
				      * and no stop() was issued in
				      * between.
				      *
				      * For some reason (error in
				      * operating system?) the
				      * function call
				      * <tt>getrusage(RUSAGE_CHILDREN,.)</tt>
				      * gives always 0 (at least
				      * on Solaris7). Hence the
				      * Timer class still does not
				      * yet work for multithreading
				      * mode.
				      */
    double              start_time_children;

				     /**
				      * Value of the wall time when start()
				      * was called the last time or when the
				      * object was created and no stop() was
				      * issued in between.
				      */
    double              start_wall_time;

				     /**
				      * Accumulated time for all previous
				      * start()/stop() cycles. The time for
				      * the present cycle is not included.
				      */
    double              cumulative_time;

				     /**
				      * Accumulated wall time for all
				      * previous start()/stop() cycles. The
				      * wall time for the present cycle is
				      * not included.
				      */
    double              cumulative_wall_time;

				     /**
				      * Store whether the timer is presently
				      * running.
				      */
    bool                running;

#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
				     /**
				      * Store whether the timer is presently
				      * running.
				      */
    MPI_Comm            mpi_communicator;

				     /**
				      * Store whether the wall time is
				      * synchronized between machines.
				      */
    bool sync_wall_time;

				     /**
				      * A structure for parallel wall time
				      * measurement that includes the minimum
				      * time recorded among all processes, the
				      * maximum time as well as the average
				      * time defined as the sum of all
				      * individual times divided by the number
				      * of MPI processes in the MPI_Comm.
				      */
    Utilities::System::MinMaxAvg mpi_data;
#endif
};



//TODO: The following class is not thread-safe
/**
 * This class can be used to generate formatted output from time
 * measurements of different subsections in a program. It is possible to
 * create several sections that perform certain aspects of the program. A
 * section can be entered several times. By changing the options in
 * OutputFrequency and OutputType, the user can choose whether output should
 * be generated every time a section is joined or just in the end of the
 * program. Moreover, it is possible to show CPU times, wall times or both.
 *
 * <h3>Usage</h3>
 *
 * Use of this class could be as follows:
 * @code
 *   TimerOuput timer (std::cout, TimerOutput::summary,
 *                     TimerOutput::wall_times);
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
 * The output will see that we entered the assembly and solve section
 * twice, and reports how much time we spent there. Moreover, the class
 * measures the total time spent from start to termination of the TimerOutput
 * object. In this case, we did a lot of other stuff, so that the time
 * proportions of the functions we measured are far away from 100 precent.
 *
 * See the step-32 tutorial program for usage of this class.
 *
 * @ingroup utilities
 * @author M. Kronbichler, 2009.
 */
class TimerOutput
{
  public:
				     /**
				      * Sets whether to generate output every
				      * time we exit a section, just in the
				      * end, or both.
				      */
    enum OutputFrequency {every_call, summary, every_call_and_summary}
    output_frequency;

				     /**
				      * Sets whether to show CPU times, wall
				      * times, or both CPU and wall times.
				      */
    enum OutputType      {cpu_times, wall_times, cpu_and_wall_times}
    output_type;

				     /**
				      * Constructor that takes std::cout as
				      * output stream.
				      */
    TimerOutput (std::ostream              &stream,
		 const enum OutputFrequency output_frequency,
		 const enum OutputType      output_type);

				     /**
				      * Constructor that takes a
				      * ConditionalOStream to write output to.
				      */
    TimerOutput (ConditionalOStream        &stream,
		 const enum OutputFrequency output_frequency,
		 const enum OutputType      output_type);

#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
				     /**
				      * Constructor that takes an MPI
				      * communicator as input. A timer
				      * constructed this way will sum up the
				      * CPU times over all processors in the
				      * MPI network for calculating the CPU
				      * time.
				      *
				      * Meant for using std::cout as output
				      * stream.
				      */
    TimerOutput (MPI_Comm                   mpi_comm,
		 std::ostream              &stream,
		 const enum OutputFrequency output_frequency,
		 const enum OutputType      output_type);

				     /**
				      * Constructor that takes an MPI
				      * communicator as input. A timer
				      * constructed this way will sum up the
				      * CPU times over all processors in the
				      * MPI network for calculating the CPU
				      * time.
				      *
				      * Constructor that takes a
				      * ConditionalOStream to write output to.
				      */
    TimerOutput (MPI_Comm                   mpi_comm,
		 ConditionalOStream        &stream,
		 const enum OutputFrequency output_frequency,
		 const enum OutputType      output_type);




#endif

				     /**
				      * Destructor. Calls print_summary() in
				      * case the option for writing the
				      * summary output is set.
				      */
    ~TimerOutput();

				     /**
				      * Open a section by given a string name
				      * of it. In case the name already
				      * exists, that section is done once
				      * again.
				      */
    void enter_subsection (const std::string &section_name);

				     /**
				      * Same as @p enter_subsection.
				      */
    void enter_section (const std::string &section_name);

				     /**
				      * Leave a section. If no name is given,
				      * the last section that was entered is
				      * left.
				      */
    void leave_subsection (const std::string &section_name = std::string());

				     /**
				      * Same as @p leave_subsection.
				      */
    void exit_section (const std::string &section_name = std::string());

				     /**
				      * Print a formatted table that
				      * summarizes the time consumed in the
				      * various sections.
				      */
    void print_summary () const;

				     /**
				      * By calling this function, all output
				      * can be disabled. This function
				      * together with enable_output() can be
				      * useful if one wants to control the
				      * output in a flexible way without
				      * putting a lot of <tt>if</tt> clauses
				      * in the program.
				      */
    void disable_output ();

				     /**
				      * This function re-enables output of
				      * this class if it was previously
				      * disabled with disable_output(). This
				      * function together with
				      * disable_output() can be useful if
				      * one wants to control the output in a
				      * flexible way without putting a lot
				      * of <tt>if</tt> clauses in the
				      * program.
				      */
    void enable_output ();

  private:
				     /**
				      * A timer object for the overall
				      * run time. If we are using MPI,
				      * this timer also accumulates
				      * over all MPI processes.
				      */
    Timer              timer_all;

				     /**
				      * A structure that groups all
				      * information that we collect
				      * about each of the sections.
				      */
    struct Section
    {
	Timer  timer;
	double total_cpu_time;
	double total_wall_time;
	unsigned int n_calls;
    };

				     /**
				      * A list of all the sections and
				      * their information.
				      */
    std::map<std::string, Section> sections;

				     /**
				      * The stream object to which we
				      * are to output.
				      */
    ConditionalOStream out_stream;

				     /**
				      * A boolean variable that sets whether
				      * output of this class is currently on
				      * or off.
				      */
    bool output_is_enabled;

				     /**
				      * A list of the sections that
				      * have been entered and not
				      * exited. The list is kept in
				      * the order in which sections
				      * have been entered, but
				      * elements may be removed in the
				      * middle if an argument is given
				      * to the exit_section()
				      * function.
				      */
    std::list<std::string> active_sections;

#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
				     /**
				      * Store whether the timer is presently
				      * running.
				      */
    MPI_Comm            mpi_communicator;
#endif

				     /**
				      * A lock that makes sure that this
				      * class gives reasonable results even
				      * when used with several threads.
				      */
    Threads::ThreadMutex mutex;
};



/* ---------------- inline functions ----------------- */


inline
void Timer::restart ()
{
  reset();
  start();
}



#ifdef DEAL_II_COMPILER_SUPPORTS_MPI

inline
const Utilities::System::MinMaxAvg &
Timer::get_data() const
{
  return mpi_data;
}



template <class STREAM>
inline
void
Timer::print_data(STREAM & stream) const
{
  unsigned int my_id = dealii::Utilities::System::get_this_mpi_process(mpi_communicator);
  if (my_id==0)
    stream << mpi_data.max << " wall,"
	   << " max @" << mpi_data.max_index
	   << ", min=" << mpi_data.min << " @" << mpi_data.min_index
	   << ", avg=" << mpi_data.avg
	   << std::endl;
}

#endif



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


DEAL_II_NAMESPACE_CLOSE

#endif
