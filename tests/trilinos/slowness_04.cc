// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2013 by the deal.II authors
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
#include <deal.II/base/utilities.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_vector.h>
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
        const unsigned int k = Testing::rand() % (N-i);
        permutation[i] = unused_indices[k];

        // then swap this used element to the
        // end where we won't consider it any
        // more
        std::swap (unused_indices[k], unused_indices[N-i-1]);
      }
  }

  // build the sparse matrix
  Epetra_Map map (TrilinosWrappers::types::int_type(N*N), 0,
                  Utilities::Trilinos::comm_world());
  TrilinosWrappers::SparseMatrix matrix (map, 5);
  for (unsigned int i_=0; i_<N; i_++)
    for (unsigned int j_=0; j_<N; j_++)
      {
        const unsigned int i=permutation[i_];
        const unsigned int j=permutation[j_];

        const unsigned int global = i*N+j;
        matrix.set(global, global, Testing::rand());
        if (j>0)
          {
            matrix.set(global-1, global, Testing::rand());
            matrix.set(global, global-1, Testing::rand());
          }
        if (j<N-1)
          {
            matrix.set(global+1, global, Testing::rand());
            matrix.set(global, global+1, Testing::rand());
          }
        if (i>0)
          {
            matrix.set(global-N, global, Testing::rand());
            matrix.set(global, global-N, Testing::rand());
          }
        if (i<N-1)
          {
            matrix.set(global+N, global, Testing::rand());
            matrix.set(global, global+N, Testing::rand());
          }
      }
  matrix.compress (VectorOperation::insert);

  // then do a single matrix-vector
  // multiplication with subsequent formation
  // of the matrix norm
  TrilinosWrappers::MPI::Vector v1(map);
  TrilinosWrappers::MPI::Vector v2(map);
  for (unsigned int i=0; i<N*N; ++i)
    v1(i) = i;
  matrix.vmult (v2, v1);

  deallog << v1 *v2 << std::endl;
}



int main (int argc,char **argv)
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv);

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
