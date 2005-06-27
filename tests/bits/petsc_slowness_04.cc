//----------------------------  petsc_slowness_03.cc  ---------------------------
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
//----------------------------  petsc_slowness_03.cc  ---------------------------


// this is part of a whole suite of tests that checks the relative speed of
// using PETSc for sparse matrices as compared to the speed of our own
// library. the tests therefore may not all actually use PETSc, but they are
// meant to compare it
//
// the tests build the 5-point stencil matrix for a uniform grid of size N*N
//
// this test does the same as the _03 test, except that it does not allocate
// the entries in consecutive order, but in pseudo-random order. the reason is
// that in usual finite element programs, we do not write to the elements of a
// matrix in a consecutive fashion, but rather according to the order of
// degrees of freedom in the sequence of cells that we traverse

#include "../tests.h"
#include <lac/sparse_matrix.h>
#include <lac/petsc_parallel_sparse_matrix.h>
#include <lac/petsc_parallel_vector.h>
#include <fstream>
#include <iostream>


void test ()
{
  const unsigned int N = 200;

                                   // first find a random permutation of the
                                   // indices
  std::vector<unsigned int> permutation (N);
  {
    std::vector<unsigned int> unused_indices (N);
    for (unsigned int i=0; i<N; i++)
      unused_indices[i] = i;

    for (unsigned int i=0; i<N; i++)
      {
                                         // pick a random element among the
                                         // unused indices
        const unsigned int k = rand() % (N-i);
        permutation[i] = unused_indices[k];
        
                                         // then swap this used element to the
                                         // end where we won't consider it any
                                         // more
        std::swap (unused_indices[k], unused_indices[N-i-1]);
      }    
  }
  
                                   // build the sparse matrix 
  PETScWrappers::MPI::SparseMatrix matrix (PETSC_COMM_WORLD,
                                           N*N, N*N,
                                           N*N, N*N,
                                           5);
  for(unsigned int i_=0; i_<N; i_++)
    for(unsigned int j_=0; j_<N; j_++)
      {
        const unsigned int i=permutation[i_];
        const unsigned int j=permutation[j_];
        
        const unsigned int global = i*N+j;
        matrix.add(global, global, rand());
        if (j>0)
          {
            matrix.add(global-1, global, rand());
            matrix.add(global, global-1, rand());
          }
        if (j<N-1)
          {
            matrix.add(global+1, global, rand());
            matrix.add(global, global+1, rand());
          }
        if (i>0)
          {
            matrix.add(global-N, global, rand());
            matrix.add(global, global-N, rand());
          }
        if (i<N-1)
          {
            matrix.add(global+N, global, rand());
            matrix.add(global, global+N, rand());
          }
      }
  matrix.compress ();
  
                                   // then do a single matrix-vector
                                   // multiplication with subsequent formation
                                   // of the matrix norm
  PETScWrappers::MPI::Vector v1(PETSC_COMM_WORLD, N*N, N*N);
  PETScWrappers::MPI::Vector v2(PETSC_COMM_WORLD, N*N, N*N);
  for (unsigned int i=0; i<N*N; ++i)
    v1(i) = i;
  matrix.vmult (v2, v1);
  
  deallog << v1*v2 << std::endl;
}



int main (int argc,char **argv) 
{
  std::ofstream logfile("petsc_slowness_04.output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  try
    {
      PetscInitialize(&argc,&argv,0,0);
      {
        test ();
      }
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
