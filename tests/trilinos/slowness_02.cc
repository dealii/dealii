//----------------------------  trilinos_slowness_02.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  trilinos_slowness_02.cc  ---------------------------


// this is part of a whole suite of tests that checks the relative speed of
// using PETSc for sparse matrices as compared to the speed of our own
// library. the tests therefore may not all actually use PETSc, but they are
// meant to compare it
//
// the tests build the 5-point stencil matrix for a uniform grid of size N*N

#include "../tests.h" 
#include <base/utilities.h>
#include <lac/sparse_matrix.h>
#include <lac/trilinos_sparse_matrix.h>
#include <lac/trilinos_vector.h>
#include <fstream>
#include <iostream>


void test ()
{
  const unsigned int N = 200;

                                   // build the sparse matrix 
  TrilinosWrappers::SparseMatrix matrix (N*N, N*N, 5U);
  for(unsigned int i=0; i<N; i++)
    for(unsigned int j=0; j<N; j++)
      {
        const unsigned int global = i*N+j;
        matrix.set(global, global, 4);
        if (j>0)
          {
            matrix.set(global-1, global, -1);
            matrix.set(global, global-1, -1);
          }
        if (j<N-1)
          {
            matrix.set(global+1, global, -1);
            matrix.set(global, global+1, -1);
          }
        if (i>0)
          {
            matrix.set(global-N, global, -1);
            matrix.set(global, global-N, -1);
          }
        if (i<N-1)
          {
            matrix.set(global+N, global, -1);
            matrix.set(global, global+N, -1);
          }
      }
  matrix.compress ();
  
  
                                   // then do a single matrix-vector
                                   // multiplication with subsequent formation
                                   // of the matrix norm
  TrilinosWrappers::Vector v1(N*N), v2(N*N);
  for (unsigned int i=0; i<N*N; ++i)
    v1(i) = i;
  matrix.vmult (v2, v1);
  
  deallog << v1*v2 << std::endl;
}



int main (int argc,char **argv) 
{
  std::ofstream logfile("slowness_02/output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  Utilities::System::MPI_InitFinalize mpi_initialization(argc, argv);

  try
    {
      {
        test ();
      }
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
  return 0;
}
