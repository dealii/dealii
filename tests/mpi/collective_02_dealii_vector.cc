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



// check Utilities::MPI::sum() for dealii::Vector

#include <deal.II/base/utilities.h>

#include <deal.II/lac/vector.h>

#include "../tests.h"

void
test()
{
  unsigned int       myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  const unsigned int numprocs = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);

  {
    Vector<int> values(2);
    values[0] = 1;
    values[1] = 2;
    Vector<int> sums(2);
    Utilities::MPI::sum(values, MPI_COMM_WORLD, sums);
    Assert((unsigned int)sums[0] == numprocs, ExcInternalError());
    Assert((unsigned int)sums[1] == 2 * numprocs, ExcInternalError());

    if (myid == 0)
      deallog << sums[0] << ' ' << sums[1] << std::endl;
  }

  {
    Vector<double> values(2);
    values[0] = 1.5;
    values[1] = 2.5;
    Vector<double> sums(2);
    Utilities::MPI::sum(values, MPI_COMM_WORLD, sums);
    Assert(sums[0] == 1.5 * numprocs, ExcInternalError());
    Assert(sums[1] == 2.5 * numprocs, ExcInternalError());

    if (myid == 0)
      deallog << sums[0] << ' ' << sums[1] << std::endl;
  }
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
