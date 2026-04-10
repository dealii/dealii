// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2004 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------



// Test distributed vector operations across multiple processes with non-local
// entries to verify correct handling and compression states.
// This tests the minimal operation to insert into a parallel distributed vector
// while not setting all non-local entries.

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


  IndexSet locally_owned(3);
  IndexSet locally_relevant(3);

  // Distribution of vector entries across processes:

  // Entry Index  |  Process Distribution  |  Owned by
  // ------------ | -----------------------| ----------
  //  0           |      0, 1              |     0
  //  1           |      0, 1              |     0
  //  2           |      0, 1              |     1


  if (my_rank == 0)
    {
      locally_owned.add_range(0, 2);
      locally_relevant.add_range(0, 3);
    }
  if (my_rank == 1)
    {
      locally_owned.add_range(2, 3);
      locally_relevant.add_range(0, 3);
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
      // Set all entries, assume all are locally active
      vector_with_nonlocal_entries[0] = 1.0;
      vector_with_nonlocal_entries[1] = 1.0;
      vector_with_nonlocal_entries[2] = 1.0;
    }
  if (my_rank == 1)
    {
      // we dont set index 0, this used to lead to
      // overwriting index 0 on rank 0.
      vector_with_nonlocal_entries[1] = 2.0;
      vector_with_nonlocal_entries[2] = 2.0;
    }

  vector_with_nonlocal_entries.compress(VectorOperation::insert);

  // Expected Results:
  // Entry Index  |  Operations    |  Result
  // ------------ | -------------- | -------
  //  0           |   =1           |   1.0
  //  1           |   =1,=1        |   1.0
  //  2           |   =2,=2        |   2.0

  // Print the result on each process
  deallog << "After operation:" << std::endl;
  vector_with_nonlocal_entries.print(deallog.get_file_stream());

  deallog << "OK" << std::endl;
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());
  MPILogInitAll log;

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
