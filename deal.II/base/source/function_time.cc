//----------------------------  function_time.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  function_time.cc  ---------------------------


#include <base/function_time.h>


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


