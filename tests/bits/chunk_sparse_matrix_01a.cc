// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2014 by the deal.II authors
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



// check setting elements in a sparse matrix using ChunkSparseMatrix::set(). make
// sure they are correct, and make sure that for the nonexisting entries
// ChunkSparseMatrix::el() returns zero and operator() throws an exception

#include "../tests.h"
#include <deal.II/lac/chunk_sparse_matrix.h>
#include <fstream>
#include <iomanip>


void test (const unsigned int chunk_size)
{
  ChunkSparsityPattern sp (5,5,3,chunk_size);
  for (unsigned int i=0; i<5; ++i)
    for (unsigned int j=0; j<5; ++j)
      if ((i+2*j+1) % 3 == 0)
        sp.add (i,j);
  sp.compress ();

  ChunkSparseMatrix<double> m(sp);

  // first set a few entries
  for (unsigned int i=0; i<m.m(); ++i)
    for (unsigned int j=0; j<m.n(); ++j)
      if ((i+2*j+1) % 3 == 0)
        m.set (i,j, i*j*.5+.5);

  // then make sure we retrieve the same ones
  for (unsigned int i=0; i<m.m(); ++i)
    for (unsigned int j=0; j<m.n(); ++j)
      if ((i+2*j+1) % 3 == 0)
        {
          Assert (m(i,j) == i*j*.5+.5, ExcInternalError());
          Assert (m.el(i,j) == i*j*.5+.5, ExcInternalError());
        }
      else
        {
          // reading elements not in the
          // sparsity pattern should return
          // zero
          const double x = m.el(i,j);
          Assert (x == 0, ExcInternalError());

          // if this is a sparsity_pattern
          // with chunk_size==1, then we need
          // to get an exception if we access
          // any other element. if
          // chunk_size>1, then this isn't
          // necessarily true. because the
          // matrix is square, we may also
          // always access the diagonal
          // element
          bool exc_thrown = false;
          try
            {
              m(i,j);
            }
          catch (const std::exception &)
            {
              exc_thrown = true;
            }
          Assert ((exc_thrown == true) || (chunk_size > 1) || (i==j),
                  ExcInternalError());
        }

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
      const unsigned int chunk_sizes[] = { 1, 2, 4, 5, 7 };
      for (unsigned int i=0;
           i<sizeof(chunk_sizes)/sizeof(chunk_sizes[0]);
           ++i)
        test (chunk_sizes[i]);
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
