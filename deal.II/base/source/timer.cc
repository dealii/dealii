/* $Id$ */

#include <basic/timer.h>
#include <ctime>
#include <sys/time.h>



#if #cpu(transputer)
#  define TIMER TimeNowLow()
#  define DIVIDE CLK_TCK_LOW
#endif
#ifdef _ARCH_PPC
#  if #system(parix)
#    define TIMER TimeNow()
#    define DIVIDE CLOCK_TICK
#  endif
#  if #system(aix)
#    define DIVIDE CLOCKS_PER_SEC
#  endif
#endif
#if #cpu(sparc)
#  define DIVIDE 1.e6
#endif

#ifndef TIMER
#define TIMER clock()
#endif

#ifndef OVER_TIME
#define OVER_TIME 4294967296./DIVIDE
#endif

#ifndef DIVIDE
#define DIVIDE 1.e6
#endif


const double Timer::overtime = OVER_TIME;


Timer::Timer()
  : cumulative_time(0.)
{
  start();
};


void Timer::start () {
  running    = true;
  overflow   = 0;
  start_time = TIMER / DIVIDE;
};



double Timer::stop () {
  running = false;
  double dtime =  TIMER / DIVIDE - start_time;
  if (dtime < 0) {
    overflow++;
  };
  
  cumulative_time += dtime;
  return full_time ();
};



double Timer::operator() () {
  if (running) {
    const double dtime =  TIMER / DIVIDE - start_time;
    if (dtime < 0) {
      overflow++;
    };
    
    return dtime + full_time();
  };
  
  return full_time();
};



void Timer::reset () {
  cumulative_time = 0.;
  running         = false;
};



double Timer::full_time () const {
  return cumulative_time + overflow*overtime;
};


