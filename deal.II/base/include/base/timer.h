//----------------------------  timer.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001 by the deal authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  timer.h  ---------------------------
#ifndef __deal2__timer_h
#define __deal2__timer_h


/**
 * This is a very simple class which provides information about the time
 * elapsed since the timer was started last time. Information is retrieved
 * from the system on the basis of clock cycles since last time the computer
 * was booted. On a SUN workstation, this information is exact to about a
 * microsecond. 
 *
 *
 * @sect3{Usage}
 *
 * Use of this class is as you might expect by looking at the member
 * functions:
 * @begin{verbatim}
 *   Time timer;
 *   timer.start ();
 *
 *   // do some complicated computations here
 *   ...
 *
 *   timer.stop ();
 *
 *   std::cout << "Elapsed time: " << timer() << " seconds.";
 *
 *   // reset timer for the next thing it shall do
 *   timer.reset();
 * @end{verbatim}
 *
 * Alternatively, you can also restart the timer instead of resetting
 * it. The times between successive calls to @p{start/stop} will then be
 * accumulated.
 *
 * Note: the implementation of this class is system dependent.
 *
 * @author G. Kanschat, W. Bangerth
 */
class Timer
{
  public:
				     /**
				      * Constructor. Starts the timer at 0 sec.
				      */
    Timer ();

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
				      * Stop the timer if neccessary and reset
				      * the elapsed time to zero.
				      */
    void reset ();

				     /**
				      * Access to the current time
				      * without disturbing time
				      * measurement. The elapsed time
				      * is returned in units of
				      * seconds.
				      */
    double operator() () const;

  private:

				     /**
				      * Value of the user time when @p{start}
				      * was called the last time or when the
				      * object was created and no @p{stop} was
				      * issued in between.
				      */
    double              start_time;


				     /**
				      * Similar to @p{start_time}, but
				      * needed for children threads
				      * in multithread mode. Value of
				      * the user time when @p{start}
				      * was called the last time or
				      * when the object was created
				      * and no @p{stop} was issued in
				      * between.
				      *
				      * For some reason (error in
				      * operating system?) the
				      * function call
				      * @p{getrusage(RUSAGE_CHILDREN,.)}
				      * gives always 0 (at least
				      * on Solaris7). Hence the
				      * @p{Timer} class still does not
				      * yet work for multithreading
				      * mode.
				      */
    double              start_time_children;

				     /**
				      * Accumulated time for all previous
				      * @p{start}/@p{stop} cycles. The time for
				      * the present cycle is not included.
				      */
    double              cumulative_time;


				     /**
				      * Store whether the timer is presently
				      * running.
				      */
    bool                running;
};


#endif
