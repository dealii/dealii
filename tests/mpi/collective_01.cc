// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check Utilities::MPI::min_max_avg()

#include <deal.II/base/utilities.h>

#include "../tests.h"

void
print_it(Utilities::MPI::MinMaxAvg &result)
{
  deallog << "sum: " << result.sum << " avg: " << result.avg
          << " min: " << result.min << " @" << result.min_index
          << " max: " << result.max << " @" << result.max_index << std::endl;
}

void
test()
{
  Assert(Utilities::MPI::job_supports_mpi(), ExcInternalError());

  unsigned int       myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  const unsigned int numprocs = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);

  if (myid == 0)
    deallog << "Running on " << numprocs << " CPU(s)." << std::endl;

  Utilities::MPI::MinMaxAvg result;

  double value = 0.0;
  result       = Utilities::MPI::min_max_avg(value, MPI_COMM_WORLD);
  if (myid == 0)
    print_it(result);
  Assert(result.sum == 0.0, ExcInternalError());
  Assert(result.min == 0.0, ExcInternalError());
  Assert(result.max == 0.0, ExcInternalError());
  Assert(result.avg == 0.0, ExcInternalError());

  value  = 1.0;
  result = Utilities::MPI::min_max_avg(value, MPI_COMM_WORLD);
  if (myid == 0)
    print_it(result);
  Assert(result.sum == numprocs, ExcInternalError());
  Assert(result.min == 1.0, ExcInternalError());
  Assert(result.max == 1.0, ExcInternalError());
  Assert(result.avg == 1.0, ExcInternalError());

  value = 0.0;
  if (myid == 0)
    value = 1.0;

  result = Utilities::MPI::min_max_avg(value, MPI_COMM_WORLD);
  if (myid == 0)
    print_it(result);
  Assert(result.sum == 1.0, ExcInternalError());
  Assert(result.min == (numprocs > 1) ? 0.0 : 1.0, ExcInternalError());
  Assert(result.max == 1.0, ExcInternalError());
  Assert(result.max_index == 0, ExcInternalError());

  if (myid == 0)
    deallog << "done" << std::endl;
}


int
main(int argc, char *argv[])
{
#ifdef DEAL_II_WITH_MPI
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());
#else
  (void)argc;
  (void)argv;
  compile_time_error;

#endif

  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    {
      initlog();

      deallog.push("mpi");
      test();
      deallog.pop();
    }
  else
    test();
}
