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


// Test Table::replicate_across_communicator(). Cycle through all
// possible root ranks.

#include <deal.II/base/table.h>

#include "../tests.h"


void
test(const unsigned int root)
{
  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == root)
    deallog << root << " of " << Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD)
            << ": ";

  Table<2, int> table(10, 12);
  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == root)
    {
      table.reinit(6, 8);
      unsigned int counter = 0;
      for (auto &i : table)
        i = counter++;
    }

  table.replicate_across_communicator(MPI_COMM_WORLD, root);

  // Now check that the sizes are correct on all processes
  Assert(table.size()[0] == 6, ExcInternalError());
  Assert(table.size()[1] == 8, ExcInternalError());

  // Then also check the table contents:
  {
    unsigned int counter = 0;

    for (auto &i : table)
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
