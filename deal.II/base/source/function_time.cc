// $Id$

#include <base/functiontime.h>


FunctionTime::FunctionTime(double initial_time)
		:
		time(initial_time)
{}

FunctionTime::~FunctionTime()
{}

void
FunctionTime::set_time (const double new_time)
{
  time = new_time;
};

void
FunctionTime::advance_time (const double delta_t)
{
  set_time (time+delta_t);
};



