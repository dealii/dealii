/* $Id$ */
/* Copyright W. Bangerth, University of Heidelberg, 1998 */


#include <base/timer.h>



// these includes should probably be properly
// ./configure'd using the AC_HEADER_TIME macro:
#include <sys/resource.h>
#include <sys/time.h>

Timer::Timer()
  : cumulative_time (0.)
{
  start();
}


void Timer::start ()
{
  running    = true;
  rusage usage;
  getrusage (RUSAGE_SELF, &usage);
  start_time = usage.ru_utime.tv_sec + 1.e-6 * usage.ru_utime.tv_usec;
}



double Timer::stop ()
{
  if (running)
    {
      running = false;
      rusage usage;
      getrusage (RUSAGE_SELF, &usage);
      const double dtime = usage.ru_utime.tv_sec + 1.e-6 * usage.ru_utime.tv_usec;
      cumulative_time += dtime - start_time;
    }
  return full_time ();
}



double Timer::operator() () const
{
  if (running)
    {
      rusage usage;
      getrusage (RUSAGE_SELF, &usage);
      const double dtime =  usage.ru_utime.tv_sec + 1.e-6 * usage.ru_utime.tv_usec;

      return dtime - start_time + full_time();
    }
  else
    return full_time();
}



void Timer::reset ()
{
  cumulative_time = 0.;
  running         = false;
}



double Timer::full_time () const
{
  return cumulative_time;
}
