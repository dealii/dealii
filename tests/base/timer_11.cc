// ---------------------------------------------------------------------
//
// Copyright (C) 2020 by the deal.II authors
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


// test that when throwing an exception inside a Timer::Scope the
// desctruction of the timer calls MPI_Abort

#include <deal.II/base/timer.h>

#include <algorithm>
#include <sstream>
#include <thread>
#include <chrono>

#include "../tests.h"

void
test()
{
  std::stringstream ss;

  TimerOutput t(MPI_COMM_WORLD, ss, TimerOutput::never, TimerOutput::wall_times);

  t.enter_subsection("Base section");

  {
  TimerOutput::Scope(t, "Test section");

  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
  {
    std::this_thread::sleep_for(std::chrono::seconds(2));
    AssertThrow(false,ExcInternalError());
  }

  MPI_Barrier(MPI_COMM_WORLD);
  }
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi(argc, argv);

  mpi_initlog();
  
  test();
}
