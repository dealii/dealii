// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------



// PETScWrappers::MPI::SparseMatrix::reinit(CompressedSparsityPattern) should
// create a matrix that, when filled with elements that match the sparsity
// pattern, doesn't require any more memory allocation any more. This is
// tricky to get right, though, and took a while until it worked.
//
// the test itself can only check this when run in parallel, but we can add it
// anyway to the testsuite. in order to check that the library actually does
// what it is supposed to do, run this program with more than one MPI process
// and with the option -log_info. All the output should say that no additional
// malloc calls have been performed

#include "../tests.h"
#include <deal.II/lac/petsc_parallel_sparse_matrix.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>
#include <fstream>
#include <cstdlib>


unsigned int
get_n_mpi_processes ()
{
  int n_jobs;
  MPI_Comm_size (MPI_COMM_WORLD, &n_jobs);

  return n_jobs;
}



unsigned int
get_this_mpi_process ()
{
  int rank;
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);

  return rank;
}


void test ()
{
  typedef PETScWrappers::MPI::SparseMatrix::size_type size_type;

  // create a parallel matrix where the first
  // process has 10 rows, the second one 20,
  // the third one 30, and so on
  unsigned int N = 0;
  std::vector<size_type> local_rows_per_process (get_n_mpi_processes());
  std::vector<size_type> start_row (get_n_mpi_processes());
  for (unsigned int i=0; i<get_n_mpi_processes(); ++i)
    {
      N += (i+1)*10;
      local_rows_per_process[i] = (i+1)*10;
      start_row[i] += i*10;
    }

  // here is a sparsity pattern for which we
  // used to allocate additional memory for 2
  // processes. note that only one of the
  // four blocks uses Inodes
  CompressedSparsityPattern csp (N,N);
  for (unsigned int i=0; i<N; ++i)
    for (unsigned int j=0; j<N; ++j)
      {
        csp.add (i,i);
        if (i+local_rows_per_process.back() < N)
          csp.add (i,i+local_rows_per_process.back());
        if (i > local_rows_per_process.back())
          csp.add (i,i-local_rows_per_process.back());
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
  PETScWrappers::MPI::SparseMatrix m;
  m.reinit (MPI_COMM_WORLD, csp, local_rows_per_process,
            local_rows_per_process, get_this_mpi_process());

  // no write into the exact same matrix
  // entries as have been created by the
  // sparsity pattern above
  for (unsigned int i=0; i<N; ++i)
    for (unsigned int j=0; j<csp.row_length(i); ++j)
      m.add (i, csp.column_number(i,j), 1.);

  m.compress (VectorOperation::add);

  deallog << "OK" << std::endl;
}



int main (int argc,char **argv)
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      test ();

    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
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
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    };
}
