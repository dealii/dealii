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



// comparisons between chunk sparse matrix iterators, same as
// sparse_matrix_iterators_11 otherwise.

#include "../tests.h"
#include <deal.II/lac/chunk_sparse_matrix.h>
#include <fstream>
#include <iomanip>


void test (const unsigned int chunk_size)
{
  deallog << "Chunk size: " << chunk_size << std::endl;

  // create a sparsity pattern with totally
  // empty lines (not even diagonals, since
  // not quadratic)
  ChunkSparsityPattern sparsity(4,5,1,chunk_size);
  sparsity.add (1,1);
  sparsity.add (3,1);
  sparsity.compress ();

  // attach a sparse matrix to it
  ChunkSparseMatrix<double> A(sparsity);

  ChunkSparseMatrix<double>::iterator k = A.begin(),
                                      j = ++A.begin();

  Assert (k < j, ExcInternalError());
  Assert (j > k, ExcInternalError());

  Assert (!(j < k), ExcInternalError());
  Assert (!(k > j), ExcInternalError());

  Assert (k != j, ExcInternalError());
  Assert (!(k == j), ExcInternalError());

  Assert (k == k, ExcInternalError());
  Assert (!(k != k), ExcInternalError());

  deallog << "OK" << std::endl;
}



int main ()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  try
    {
      test (1);
      test (3);
    }
  catch (std::exception &exc)
    {
      deallog << std::endl << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Exception on processing: " << std::endl
              << exc.what() << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;

      return 1;
    }
  catch (...)
    {
      deallog << std::endl << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Unknown exception!" << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
      return 1;
    };
}
