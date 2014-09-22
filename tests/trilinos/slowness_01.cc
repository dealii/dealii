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
// using Trilinos for sparse matrices as compared to the speed of our own
// library. the tests therefore may not all actually use Trilinos, but they are
// meant to compare it
//
// the tests build the 5-point stencil matrix for a uniform grid of size N*N

#include "../tests.h"
#include <deal.II/base/utilities.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <fstream>
#include <iostream>


void test ()
{
  const unsigned int N = 200;

  // first build the sparsity pattern
  SparsityPattern sparsity (N*N, N*N, 5);
  for (unsigned int i=0; i<N; i++)
    for (unsigned int j=0; j<N; j++)
      {
        const unsigned int global = i*N+j;
        sparsity.add(global, global);
        if (j>0)
          {
            sparsity.add(global-1, global);
            sparsity.add(global, global-1);
          }
        if (j<N-1)
          {
            sparsity.add(global+1, global);
            sparsity.add(global, global+1);
          }
        if (i>0)
          {
            sparsity.add(global-N, global);
            sparsity.add(global, global-N);
          }
        if (i<N-1)
          {
            sparsity.add(global+N, global);
            sparsity.add(global, global+N);
          }
      }
  sparsity.compress();


  // next build the sparse matrix itself
  SparseMatrix<double> matrix (sparsity);
  for (unsigned int i=0; i<N; i++)
    for (unsigned int j=0; j<N; j++)
      {
        const unsigned int global = i*N+j;
        matrix.add(global, global, 4);
        if (j>0)
          {
            matrix.add(global-1, global, -1);
            matrix.add(global, global-1, -1);
          }
        if (j<N-1)
          {
            matrix.add(global+1, global, -1);
            matrix.add(global, global+1, -1);
          }
        if (i>0)
          {
            matrix.add(global-N, global, -1);
            matrix.add(global, global-N, -1);
          }
        if (i<N-1)
          {
            matrix.add(global+N, global, -1);
            matrix.add(global, global+N, -1);
          }
      }

  // then do a single matrix-vector
  // multiplication with subsequent formation
  // of the matrix norm
  Vector<double> v1(N*N), v2(N*N);
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

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv);

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
