// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check if mpi is working

#include <deal.II/base/mpi.h>
#include <deal.II/base/utilities.h>

#include "../tests.h"


void
test_mpi()
{
  Assert(Utilities::MPI::job_supports_mpi(), ExcInternalError());


  unsigned int       myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  const unsigned int numprocs = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);

  if (myid == 0)
    deallog << "Running on " << numprocs << " CPU(s)." << std::endl;

  for (unsigned int i = 1; i < numprocs; ++i)
    {
      MPI_Barrier(MPI_COMM_WORLD);

      if (myid == 0)
        {
          unsigned int buf = numbers::invalid_unsigned_int;
          MPI_Status   status;
          MPI_Recv(&buf, 1, MPI_UNSIGNED, i, 1, MPI_COMM_WORLD, &status);
          deallog << "got message '" << buf << "' from CPU " << i + 1 << '!'
                  << std::endl;
          Assert(buf == i, ExcInternalError());
        }
      else if (myid == i)
        {
          MPI_Send(&myid, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
        }
    }
  if (myid == 0)
    deallog << "done" << std::endl;
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());

  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    {
      initlog();

      deallog.push("mpi");
      test_mpi();
      deallog.pop();
    }
  else
    test_mpi();
}
