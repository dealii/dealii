// ---------------------------------------------------------------------
//
// Copyright (C) 2012 - 2018 by the deal.II authors
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
