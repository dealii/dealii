// ---------------------------------------------------------------------
//
// Copyright (C) 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

// Basic tests for class TimestepControl

#include <deal.II/algorithms/timestep_control.h>

#include "../tests.h"

using Algorithms::TimestepControl;

void
print_time(TimestepControl &time)
{
  deallog << "Current time = " << time.now() << ", step size = " << time.step()
          << ", Generate output = " << (time.print() ? "yes" : "no")
          << std::endl;
}

void
test_uniform_strategy()
{
  LogStream::Prefix p("uniform");

  TimestepControl time(/*start*/ 0.,
                       /*final*/ 1.5,
                       /*tolerance*/ 0.1,
                       /*start_step*/ 0.123,
                       /*print_step*/ 0.31);
  time.strategy(TimestepControl::uniform);

  print_time(time);
  while (time.now() != time.final())
    {
      time.advance();
      print_time(time);
    }
}

void
test_doubling_strategy()
{
  LogStream::Prefix p("doubling");

  TimestepControl time;
  // using the set methods to configure the object
  time.start(0.);
  time.final(2.53);
  time.tolerance(0.1);
  time.start_step(0.123);
  time.strategy(TimestepControl::doubling);
  time.restart();

  print_time(time);
  while (time.now() != time.final())
    {
      time.advance();
      print_time(time);
    }
}

int
main()
{
  initlog();
  test_uniform_strategy();
  test_doubling_strategy();
  deallog << "OK" << std::endl;
  return 0;
}
