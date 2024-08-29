// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2004 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Test distributed vector operations across multiple processes with non-local
// entries to verify correct handling and compression states.

#include <deal.II/base/index_set.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/trilinos_tpetra_vector.h>

#include <iostream>
#include <vector>

#include "../tests.h"


void
test()
{
  unsigned int my_rank = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  unsigned int n_procs = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);


  IndexSet locally_owned(8);
  IndexSet locally_relevant(8);

  // Distribution of vector entries across processes:
  // Process 0 and Process 1 have non-local entries, while Process 2 does not.

  // Entry Index  |  Process Distribution  |  Owned by
  // ------------ | ----------------------| ----------
  //  0           |      0                 |     0
  //  1           |      0                 |     0
  //  2           |      0, 1              |     0
  //  3           |      0, 1              |     1
  //  4           |         1              |     1
  //  5           |         1, 2           |     2
  //  6           |            2           |     2
  //  7           |            2           |     2


  if (my_rank == 0)
    {
      locally_owned.add_range(0, 3);
      locally_relevant.add_range(0, 4);
    }
  if (my_rank == 1)
    {
      locally_owned.add_range(3, 5);
      locally_relevant.add_range(2, 6);
    }
  if (my_rank == 2)
    {
      // no nonlocal entries on rank 2
      locally_owned.add_range(5, 8);
      locally_relevant.add_range(5, 8);
    }
  locally_owned.compress();
  locally_relevant.compress();

  // Create a vector with writable entries
  LinearAlgebra::TpetraWrappers::Vector<double, MemorySpace::Default>
    vector_with_nonlocal_entries;

  vector_with_nonlocal_entries.reinit(locally_owned,
                                      locally_relevant,
                                      MPI_COMM_WORLD,
                                      true);

  vector_with_nonlocal_entries = 0.0;

  if (my_rank == 0)
    {
      // Add to both local and non-local parts
      vector_with_nonlocal_entries[0] += 1.0;
      vector_with_nonlocal_entries[1] += 1.0;
      vector_with_nonlocal_entries[2] += 1.0;
      vector_with_nonlocal_entries[3] += 1.0;
    }
  if (my_rank == 1)
    {
      // Add to both local and non-local parts
      vector_with_nonlocal_entries[2] += 2.0;
      vector_with_nonlocal_entries[3] += 2.0;
      vector_with_nonlocal_entries[4] += 2.0;
      vector_with_nonlocal_entries[5] += 2.0;
    }
  if (my_rank == 2)
    {
      // Add to only local part, there's no nonlocal part on rank 2
      vector_with_nonlocal_entries[5] += 3.0;
      vector_with_nonlocal_entries[6] += 3.0;
      vector_with_nonlocal_entries[7] += 3.0;
    }
  vector_with_nonlocal_entries.compress(VectorOperation::add);

  // Expected Results:
  // Entry Index  |  Additions     |  Result
  // ------------ | -------------- | -------
  //  0           |   +1           |   1.0
  //  1           |   +1           |   1.0
  //  2           |   +1, +2       |   3.0
  //  3           |   +1, +2       |   3.0
  //  4           |       +2       |   2.0
  //  5           |       +2, +3   |   5.0
  //  6           |           +3   |   3.0
  //  7           |           +3   |   3.0

  // Check the results on each process
  if (my_rank == 0)
    {
      AssertThrow((vector_with_nonlocal_entries[0] == 1.0), ExcInternalError());
      AssertThrow((vector_with_nonlocal_entries[1] == 1.0), ExcInternalError());
      AssertThrow((vector_with_nonlocal_entries[2] == 3.0), ExcInternalError());
    }
  if (my_rank == 1)
    {
      AssertThrow((vector_with_nonlocal_entries[3] == 3.0), ExcInternalError());
      AssertThrow((vector_with_nonlocal_entries[4] == 2.0), ExcInternalError());
    }
  if (my_rank == 2)
    {
      AssertThrow((vector_with_nonlocal_entries[5] == 5.0), ExcInternalError());
      AssertThrow((vector_with_nonlocal_entries[6] == 3.0), ExcInternalError());
      AssertThrow((vector_with_nonlocal_entries[7] == 3.0), ExcInternalError());
    }

  deallog << "OK" << std::endl;
}



int
main(int argc, char **argv)
{
  initlog();

  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());


  try
    {
      test();
    }
  catch (const std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;

      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    };
}
