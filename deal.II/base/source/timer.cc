/* $Id$ */
/* Copyright W. Bangerth, University of Heidelberg, 1998 */


#include <basic/timer.h>
#include <ctime>
#include <sys/time.h>



// maybe use times() instead of clock()?


const double overtime = 4294967296./CLOCKS_PER_SEC;



Timer::Timer()
  : cumulative_time(0.)
{
  start();
};


void Timer::start () {
  running    = true;
  overflow   = 0;
  start_time = static_cast<double>(clock()) /
	       CLOCKS_PER_SEC;
};



double Timer::stop () {
  running = false;
  double dtime =  (static_cast<double>(clock()) / CLOCKS_PER_SEC -
		   start_time);
  if (dtime < 0) {
    overflow++;
  };
  
  cumulative_time += dtime;
  return full_time ();
};



double Timer::operator() () {
  if (running)
    {
      const double dtime =  static_cast<double>(clock()) / CLOCKS_PER_SEC - start_time;
      if (dtime < 0)
	overflow++;

      return dtime + full_time();
    }
  else
    return full_time();
};



void Timer::reset () {
  cumulative_time = 0.;
  running         = false;
};



double Timer::full_time () const {
  return cumulative_time + overflow*overtime;
};


