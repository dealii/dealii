// ---------------------------------------------------------------------
//
// Copyright (C) 2019 - 2020 by the deal.II authors
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

// Basic tests for class DiscreteTime

#include <deal.II/base/discrete_time.h>

#include "../tests.h"

void
print_time(const DiscreteTime &time)
{
  if (time.is_at_start())
    deallog << "Simulation started" << std::endl;
  deallog << "Current time = " << time.get_current_time()
          << ", next = " << time.get_next_time()
          << ", previous = " << time.get_previous_time()
          << ", step number = " << time.get_step_number()
          << ", next step size = " << time.get_next_step_size()
          << ", previous step size = " << time.get_previous_step_size()
          << std::endl;
  if (time.is_at_end())
    deallog << "Simulation ended" << std::endl;
}

void
test_from_start_to_end()
{
  LogStream::Prefix p("Start to end");

  DiscreteTime time(/*start_time*/ 0.,
                    /*end_time*/ 1.5,
                    /*desired_start_step_size*/ 0.123);

  deallog << "Start time = " << time.get_start_time() << std::endl;
  deallog << "End time = " << time.get_end_time() << std::endl;
  print_time(time);
  while (time.get_current_time() != time.get_end_time())
    {
      time.advance_time();
      print_time(time);
    }
  time.restart();
  deallog << "Restarted" << std::endl;
  print_time(time);
  deallog << "OK" << std::endl;
}

void
test_adjust_time_step_size()
{
  LogStream::Prefix p("Adjust time step size");

  DiscreteTime time(/*start_time*/ 0.4,
                    /*end_time*/ 2.1,
                    /*desired_start_step_size*/ 0.15);
  print_time(time);
  time.advance_time();
  print_time(time);
  time.set_desired_next_step_size(0.36);
  time.advance_time();
  print_time(time);
  time.advance_time();
  print_time(time);
  time.set_desired_next_step_size(0.61);
  time.advance_time();
  print_time(time);
  time.advance_time(); // Here we reach the end time.
  print_time(time);
  // If we call time.advance_time() one more time, it fails an assersion
  deallog << "OK" << std::endl;
}

int
main()
{
  initlog();
  test_from_start_to_end();
  test_adjust_time_step_size();
  return 0;
}
