// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test TimerOutput::print_wall_time_statistics() with optional
// quantile argument

#include <deal.II/base/timer.h>

#include <algorithm>
#include <sstream>

#include "../tests.h"

// burn computer time
double s = 0.;
void
burn(unsigned int n)
{
  for (unsigned int i = 0; i < n; ++i)
    {
      for (unsigned int j = 1; j < 100000; ++j)
        {
          s += 1. / j * i;
        }
    }
}

void
test()
{
  std::stringstream ss;

  TimerOutput t(ss, TimerOutput::never, TimerOutput::wall_times);

  t.enter_subsection("hi");
  burn(50);
  t.leave_subsection("hi");

  t.print_wall_time_statistics(MPI_COMM_WORLD, 0.1);

  t.enter_subsection(
    "this is a very long section name that previously did not work");
  burn(50);
  t.leave_subsection(
    "this is a very long section name that previously did not work");

  t.print_wall_time_statistics(MPI_COMM_WORLD, 0.1);

  std::string s = ss.str();
  std::replace_if(s.begin(), s.end(), ::isdigit, ' ');
  std::replace_if(
    s.begin(), s.end(), [](char x) { return x == '.'; }, ' ');

  deallog << s << std::endl << std::endl;
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi(argc, argv);

  mpi_initlog();

  test();
}
