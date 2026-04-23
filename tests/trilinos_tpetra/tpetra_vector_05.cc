// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2018 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


#include <deal.II/base/utilities.h>

#include <deal.II/lac/read_write_vector.h>
#include <deal.II/lac/trilinos_tpetra_vector.h>

#include <iostream>
#include <vector>

#include "../tests.h"

// Check LinearAlgebra::TpetraWrappers::Vector sadd with ghosted vector

template <typename Number>
void
test()
{
  IndexSet parallel_partitioner_distributed(10);
  IndexSet parallel_partitioner_ghosted(10);

  unsigned int rank = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  if (rank == 0)
    {
      parallel_partitioner_distributed.add_range(0, 5);
      parallel_partitioner_ghosted.add_range(0, 7);
    }
  else
    {
      parallel_partitioner_distributed.add_range(5, 10);
      parallel_partitioner_ghosted.add_range(3, 10);
    }

  parallel_partitioner_distributed.compress();
  parallel_partitioner_ghosted.compress();

  LinearAlgebra::TpetraWrappers::Vector<Number, MemorySpace::Default>
    vector_distributed(parallel_partitioner_distributed, MPI_COMM_WORLD);
  LinearAlgebra::TpetraWrappers::Vector<Number, MemorySpace::Default>
    vector_ghosted(parallel_partitioner_distributed,
                   parallel_partitioner_ghosted,
                   MPI_COMM_WORLD);

  if (rank == 0)
    {
      for (unsigned int i = 0; i < 5; ++i)
        {
          vector_distributed[i] = i;
        }
    }
  else
    {
      for (unsigned int i = 5; i < 10; ++i)
        {
          vector_distributed[i] = i;
        }
    }

  vector_ghosted = vector_distributed;

  vector_distributed.sadd(1., 1., vector_ghosted);

  if (rank == 0)
    {
      for (unsigned int i = 0; i < 5; ++i)
        AssertThrow(vector_distributed[i] == 2 * i,
                    ExcMessage("Problem in sadd()."));
    }
  else
    {
      for (unsigned int i = 5; i < 10; ++i)
        AssertThrow(vector_distributed[i] == 2 * i,
                    ExcMessage("Problem in sadd()."));
    }
}


int
main(int argc, char **argv)
{
  initlog();
  deallog.depth_console(0);

  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv, 1);

  test<double>();

  deallog << "OK" << std::endl;

  return 0;
}
