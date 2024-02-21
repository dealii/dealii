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


// Test AlignedVector::replicate_across_communicator()

#include <deal.II/base/aligned_vector.h>

#include "../tests.h"


void
test()
{
  const MPI_Comm     communicator = MPI_COMM_WORLD;
  const unsigned int root         = 1;
  Assert(root < Utilities::MPI::n_mpi_processes(communicator),
         ExcInternalError());

  // Create a vector of differing sizes on all processes, but only on
  // the root process do we put something useful into it
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

  // Final step, let every process output what it now has
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
