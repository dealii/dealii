/*----------------------------   timer.h     ---------------------------*/
/*      $Id$                 */
/*      Copyright W. Bangerth, University of Heidelberg, 1998 */
#ifndef __timer_H
#define __timer_H
/*----------------------------   timer.h     ---------------------------*/


/**
 * This is a very simple class which provides information about the time
 * elapsed since the timer was started last time. Information is retrieved
 * from the system on the basis of clock cycles since last time the computer
 * was booted. On a SUN workstation, this information is exact to about a
 * microsecond. 
 *
 *
 * \subsection{Usage}
 *
 * Use of this class is as you might expect by looking at the member
 * functions:
 * \begin{verbatim}
 *   Time timer;
 *   timer.start ();
 *
 *   // do some complicated computations here
 *   ...
 *
 *   timer.stop ();
 *
 *   cout << "Elapsed time: " << timer() << " seconds.";
 *
 *   // reset timer for the next thing it shall do
 *   timer.reset();
 * \end{verbatim}
 *
 * Alternatively, you can also restart the timer instead of resetting
 * it. The times between successive calls to #start/stop# will then be
 * accumulated.
 *
 * Note: the implementation of this class is system dependant.
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
				      *
				      * A regular call to this
				      * function serves to avoid time
				      * overflow (which is now nearly
				      * every 30 minutes on UNIX
				      * machines) in long-time
				      * measurements.
				      */
    double operator() () const;

  private:

				     /**
				      * Value of the user time when #start#
				      * was called the last time or when the
				      * object was created and no #stop# was
				      * issued in between.
				      */
    double              start_time;

				     /**
				      * Accumulated time for all previous
				      * #start#/#stop# cycles. The time for
				      * the present cycle is not included.
				      */
    double              cumulative_time;


				     /**
				      * Store whether the timer is presently
				      * running.
				      */
    bool                running;

				     /**
				      * Return #cumulative_time# plus the
				      * number of overflows times the number
				      * of seconds per overflow. Do not include
				      * the time since the last overflow
				      * occured.
				      */
    double full_time () const;
};





/*----------------------------   timer.h     ---------------------------*/
/* end of #ifndef __timer_H */
#endif
/*----------------------------   timer.h     ---------------------------*/
