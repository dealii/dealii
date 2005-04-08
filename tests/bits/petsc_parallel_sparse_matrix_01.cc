//----------------------------  petsc_parallel_sparse_matrix_01.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2004, 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  petsc_parallel_sparse_matrix_01.cc  ---------------------------


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
#include <lac/petsc_parallel_sparse_matrix.h>
#include <lac/compressed_sparsity_pattern.h>
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
                                   // create a parallel matrix where the first
                                   // process has 10 rows, the second one 20,
                                   // the third one 30, and so on
  unsigned int N = 0;
  std::vector<unsigned int> local_rows_per_process (get_n_mpi_processes());
  std::vector<unsigned int> start_row (get_n_mpi_processes());
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

  m.compress ();

  deallog << "OK" << std::endl;
}



int main (int argc,char **argv) 
{
  std::ofstream logfile("petsc_parallel_sparse_matrix_01.output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  try
    {
      PetscInitialize(&argc,&argv,0,0);
      test ();
      PetscFinalize();
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
