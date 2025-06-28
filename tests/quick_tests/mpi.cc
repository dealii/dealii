// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// test that MPI is working correctly. Note that this test expects to
// be executed with exactly two threads.

#include <deal.II/base/mpi.h>

#include <deal.II/grid/tria.h>

#include <sched.h>

#include <iostream>

int
main(int argc, char *argv[])
{
  MPI_Init(&argc, &argv);

  int myrank, nproc;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  std::cout << " Hi from " << myrank << '/' << nproc << std::endl;

  if (nproc != 2)
    {
      std::cerr << "ERROR: process does not see nproc=2!" << std::endl;
      return -1;
    }

  MPI_Barrier(MPI_COMM_WORLD);

  int err   = MPI_SUCCESS;
  int value = myrank;

  if (myrank == 1)
    err = MPI_Send(&value, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
  else if (myrank == 0)
    err = MPI_Recv(&value, 1, MPI_INT, 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  if (err != MPI_SUCCESS)
    {
      std::cerr << "MPI_Send/Recv error code = " << err << std::endl;
      abort();
    }

  if (myrank == 0 && value != 1)
    {
      std::cerr << "ERROR: MPI_Send/Recv did not work!" << std::endl;
      return -1;
    }

  value      = 1;
  int output = 0;

  MPI_Allreduce(&value, &output, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if (output != nproc)
    {
      std::cerr << "ERROR: MPI_Allreduce doesn't seem to work!" << std::endl;
      return -1;
    }

  // we need this, otherwise gcc will not link against deal.II
  dealii::Triangulation<2> test;

  MPI_Finalize();
}
