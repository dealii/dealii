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



// check Utilities::MPI::min()

#include <deal.II/base/utilities.h>

#include "../tests.h"

void
test()
{
  unsigned int       myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  const unsigned int numprocs = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);

  int          int_sum;
  unsigned int uint_sum;
  double       double_sum;
  float        float_sum;

  int_sum  = Utilities::MPI::min<int>(numprocs + myid, MPI_COMM_WORLD);
  uint_sum = Utilities::MPI::min<unsigned int>(numprocs + myid, MPI_COMM_WORLD);
  float_sum  = Utilities::MPI::min<float>(numprocs + myid, MPI_COMM_WORLD);
  double_sum = Utilities::MPI::min<double>(numprocs + myid, MPI_COMM_WORLD);

  if (myid == 0)
    deallog << int_sum << ' ' << uint_sum << ' ' << double_sum << ' '
            << float_sum << std::endl;
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
