// ---------------------------------------------------------------------
// $Id$
//
// Copyright (C) 1999 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#include <deal.II/base/function_time.h>

DEAL_II_NAMESPACE_OPEN


FunctionTime::FunctionTime(const double initial_time)
  :
  time(initial_time)
{}



FunctionTime::~FunctionTime()
{}



void
FunctionTime::set_time (const double new_time)
{
  time = new_time;
}



void
FunctionTime::advance_time (const double delta_t)
{
  set_time (time+delta_t);
}

DEAL_II_NAMESPACE_CLOSE
