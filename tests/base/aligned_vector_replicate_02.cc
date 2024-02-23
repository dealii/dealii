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


// Test AlignedVector::replicate_across_communicator(). This function
// puts data into a shared memory area for all MPI processes that are
// on the same machine. This is a read-write area, and so writing into
// it becomes visible also to those processes on the same
// machine. Since we run the test suite only in setups where *all*
// processes are on the same machine, we know that a change made on
// one process should be visible on *all*.

#include <deal.II/base/aligned_vector.h>

#include "../tests.h"


void
test()
{
  const MPI_Comm     communicator = MPI_COMM_WORLD;
  const unsigned int root         = 1;
  Assert(root < Utilities::MPI::n_mpi_processes(communicator),
         ExcInternalError());

  // Same as for the _01 test:
  AlignedVector<int> avec(
    Utilities::MPI::this_mpi_process(communicator) == root ? 10 : 5);
  if (Utilities::MPI::this_mpi_process(communicator) == root)
    {
      avec[2] = 2;
      avec[4] = 4;
      avec[6] = 6;
    }
  else
    {
      for (unsigned int i = 0; i < avec.size(); ++i)
        avec[i] = 100 * Utilities::MPI::this_mpi_process(communicator) + i;
    }

  // Now replicate this information across all processes
  avec.replicate_across_communicator(communicator, root);


  // Now let each process change one element:
  if (Utilities::MPI::this_mpi_process(communicator) == 0)
    avec[0] = 42;
  else if (Utilities::MPI::this_mpi_process(communicator) == 1)
    avec[1] = 43;
  else if (Utilities::MPI::this_mpi_process(communicator) == 2)
    avec[3] = 44;

  // Create a barrier before we check the effect
  MPI_Barrier(communicator);

  // Final step, let every process output what it now has. Every
  // process should be able to see the changed element, even though
  // only one changed it.
  deallog << "On process " << Utilities::MPI::this_mpi_process(communicator)
          << ": " << std::endl;
  for (const auto i : avec)
    deallog << i << ' ';
  deallog << std::endl;
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    all;

  test();
}
