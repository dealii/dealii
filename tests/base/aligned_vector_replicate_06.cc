// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test AlignedVector::replicate_across_communicator(). This variation
// of the _03 test checks that we can move a replicated object and
// that the moved-to object inherits the shared memory window and all
// other information. This tests uses move-assignment.

#include <deal.II/base/aligned_vector.h>

#include "../tests.h"


void
test(const unsigned int root)
{
  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == root)
    deallog << root << " of " << Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD)
            << ": ";

  AlignedVector<int> avec(1000);
  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == root)
    {
      unsigned int counter = 0;
      for (auto &i : avec)
        i = counter++;
    }

  avec.replicate_across_communicator(MPI_COMM_WORLD, root);

  // Move avec
  AlignedVector<int> xvec;
  xvec = std::move(avec);
  {
    unsigned int counter = 0;

    for (auto &i : xvec)
      {
        AssertDimension(i, counter);
        ++counter;
      }
  }

  MPI_Barrier(MPI_COMM_WORLD);

  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == root)
    deallog << "OK!" << std::endl;
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    all;

  for (unsigned int i = 0; i < Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
       ++i)
    test(i);
}
