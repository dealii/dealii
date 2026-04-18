// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2004 - 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------



// Tests the PETScWrappers::MPI::FullMatrix class

#include <deal.II/lac/petsc_full_matrix.h>

#include "../tests.h"



unsigned int
get_n_mpi_processes()
{
  int n_jobs;
  MPI_Comm_size(MPI_COMM_WORLD, &n_jobs);

  return n_jobs;
}



unsigned int
get_this_mpi_process()
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  return rank;
}


void
test()
{
  using size_type = PETScWrappers::MPI::FullMatrix::size_type;

  // create a parallel matrix where the first
  // process has 10 rows, the second one 20,
  // the third one 30, and so on
  unsigned int           N = 0;
  std::vector<size_type> local_rows_per_process(get_n_mpi_processes());
  std::vector<size_type> start_row(get_n_mpi_processes());
  for (unsigned int i = 0; i < get_n_mpi_processes(); ++i)
    {
      N += (i + 1) * 10;
      local_rows_per_process[i] = (i + 1) * 10;
      start_row[i] += i * 10;
    }


  // now create a matrix
  PETScWrappers::MPI::FullMatrix m;
  m.reinit(MPI_COMM_WORLD,
           N,
           N,
           local_rows_per_process,
           local_rows_per_process,
           get_this_mpi_process());

  for (unsigned int i = 0; i < N; ++i)
    for (unsigned int j = 0; j < N; ++j)
      m.add(i, j, 1.);

  deallog << "OK" << std::endl;
}



int
main(int argc, char **argv)
{
  initlog();

  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
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
