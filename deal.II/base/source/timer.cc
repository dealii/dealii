/* $Id$ */

#include <basic/timer.h>
#include <ctime>
#include <sys/time.h>


#ifndef TIMER
#define TIMER static_cast<double>(clock())
#endif

#ifndef DIVIDE
#define DIVIDE CLOCKS_PER_SEC
#endif
#ifndef OVER_TIME
#define OVER_TIME 4294967296./DIVIDE
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


