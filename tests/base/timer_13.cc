// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2020 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


// test TimerOutput::print_wall_time_statistics() if not all ranks visit all
// subsections

#include <deal.II/base/timer.h>

#include <algorithm>
#include <chrono>
#include <sstream>
#include <thread>

#include "../tests.h"

void
test()
{
  std::stringstream ss;

  TimerOutput t(ss, TimerOutput::never, TimerOutput::wall_times);

  const int rank = dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  if (rank == 0)
    {
      t.enter_subsection("hi");
      std::this_thread::sleep_for(std::chrono::milliseconds(1));
      t.leave_subsection("hi");
    }
  // Make sure the minimal time reported does not contain 'e-'
  std::this_thread::sleep_for(std::chrono::milliseconds(1));
  t.print_wall_time_statistics(MPI_COMM_WORLD);

  t.enter_subsection(
    "this is a very long section name that previously did not work");
  std::this_thread::sleep_for(std::chrono::milliseconds(1));
  t.leave_subsection(
    "this is a very long section name that previously did not work");

  t.print_wall_time_statistics(MPI_COMM_WORLD);

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
