/*----------------------------   timer.h     ---------------------------*/
/*      $Id$                 */
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
 *
 * Note: the implementation of this class is system dependant.
 *
 * @author R. Becker, G. Kanschat, F.-T. Suttmeier, revised by W. Bangerth
 */
class Timer {
public:
				     /**
				      * Constructor. Starts the timer at 0 sec.
				      */
    Timer();

				     /**
				      * Re-start the timer at the point where
				      * it was stopped. This way a cumulative
				      * measurement of time is possible.
				      */
    void start();

				     /**
				      * Sets the current time as next starting
				      * time and return it.
				      */
    double stop();

				     /**
				      * Stop the timer if neccessary and reset
				      * the elapsed time to zero.
				      */
    void reset();

				     /**
				      * Access to the current time without
				      * disturbing time measurement.
				      */
  // A regular call to this function serves to avoid time overflow
  // (which is now nearly every 30 minutes on UNIX machines) in long-time
  // measurements.
  //
    double operator() ();

  private:
    
    double              start_time;
    double              cumulative_time;
    static const double overtime;
    unsigned int        overflow;

				     /**
				      * Store whether the timer is presently
				      * running.
				      */
    bool                running;
    
    double full_time() const;
};





/*----------------------------   timer.h     ---------------------------*/
/* end of #ifndef __timer_H */
#endif
/*----------------------------   timer.h     ---------------------------*/
