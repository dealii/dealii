// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2004 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// TrilinosWrappers::SparseMatrix::reinit(DynamicSparsityPattern) should
// create a matrix that, when filled with elements that match the sparsity
// pattern, doesn't require any more memory allocation any more. This is
// tricky to get right, though, and took a while until it worked.
//
// the test itself can only check this when run in parallel, but we can add it
// anyway to the testsuite. in order to check that the library actually does
// what it is supposed to do, run this program with more than one MPI process
// and with the option -log_info. All the output should say that no additional
// malloc calls have been performed

#include <deal.II/base/mpi.h>
#include <deal.II/base/trilinos_utilities.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>

#include "../tests.h"


void
test()
{
  // create a parallel matrix where the first
  // process has 10 rows, the second one 20,
  // the third one 30, and so on
  unsigned int              N = 0;
  std::vector<unsigned int> local_rows_per_process(
    Utilities::Trilinos::get_n_mpi_processes(
      Utilities::Trilinos::comm_world()));
  std::vector<unsigned int> start_row(Utilities::Trilinos::get_n_mpi_processes(
    Utilities::Trilinos::comm_world()));
  for (unsigned int i = 0; i < Utilities::Trilinos::get_n_mpi_processes(
                                 Utilities::Trilinos::comm_world());
       ++i)
    {
      N += (i + 1) * 10;
      local_rows_per_process[i] = (i + 1) * 10;
      start_row[i] += i * 10;
    }

  // here is a sparsity pattern for which we
  // used to allocate additional memory for 2
  // processes. note that only one of the
  // four blocks uses Inodes
  DynamicSparsityPattern csp(N, N);
  for (unsigned int i = 0; i < N; ++i)
    for (unsigned int j = 0; j < N; ++j)
      {
        csp.add(i, i);
        if (i + local_rows_per_process.back() < N)
          csp.add(i, i + local_rows_per_process.back());
        if (i > local_rows_per_process.back())
          csp.add(i, i - local_rows_per_process.back());
      }

  // here is a sparsity pattern for which no
  // Inodes are used, but it doesn't allocate
  // additional memory
  //   for (unsigned int bi=0; bi<get_n_mpi_processes(); ++bi)
  //     for (unsigned int bj=0; bj<get_n_mpi_processes(); ++bj)
  //       for (unsigned int i=0; i<local_rows_per_process[bi]; ++i)
  //         for (unsigned int k=0; k<6; ++k)
  //           csp.add (start_row[bi] + i,
  //                    start_row[bj] + (i+2*k) % local_rows_per_process[bj]);


  // now create a matrix with this sparsity
  // pattern
  const unsigned int my_id = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  IndexSet           locally_owned_dofs(N);
  locally_owned_dofs.add_range(start_row[my_id], local_rows_per_process[my_id]);

  TrilinosWrappers::SparseMatrix m;
  m.reinit(locally_owned_dofs, locally_owned_dofs, csp, MPI_COMM_WORLD);
  // now write into the exact same matrix
  // entries as have been created by the
  // sparsity pattern above
  for (unsigned int i = 0; i < N; ++i)
    for (unsigned int j = 0; j < csp.row_length(i); ++j)
      m.add(i, csp.column_number(i, j), 1.);

  m.compress(VectorOperation::add);

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
