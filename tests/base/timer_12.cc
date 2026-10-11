// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2020 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


// Test that we can safely propagate an exception that is only thrown
// on some (not all) of the MPI ranks up the call-stack without
// triggering a deadlock in the Timer, TimerOutput, or TimerOutput::Scope
// classes.

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
  TimerOutput       t(MPI_COMM_WORLD,
                ss,
                TimerOutput::never,
                TimerOutput::wall_times);
  {
    TimerOutput::Scope timer_section(t, "Test section");

    if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
      {
        std::this_thread::sleep_for(std::chrono::seconds(2));
        throw std::exception();
      }

    MPI_Barrier(MPI_COMM_WORLD);
  }
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi(argc, argv);

  mpi_initlog();

  try
    {
      deallog << "Starting test ..." << std::endl;
      test();
    }
  catch (...)
    {
      deallog
        << "SUCCESS: This program propagated an exception out of the timer scope."
        << std::endl;

      throw;
    }
}
