// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

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
  time.set_next_step_size(0.36);
  time.advance_time();
  print_time(time);
  time.advance_time();
  print_time(time);
  time.set_desired_next_step_size(0.61);
  time.advance_time();
  print_time(time);
  time.advance_time(); // Here we reach the end time.
  print_time(time);
  // If we call time.advance_time() one more time, it fails an assertion
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
