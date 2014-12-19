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



// set a few elements in a chunk sparse matrix, perform a matrix-vector
// product through the iterator and compare with a matrix-vector product

#include "../tests.h"
#include <deal.II/lac/chunk_sparse_matrix.h>
#include <deal.II/lac/vector.h>
#include <fstream>
#include <iomanip>


void test (const unsigned int chunk_size)
{
  deallog << "Chunk size = " << chunk_size << std::endl;

  for (unsigned int n_cols = 4; n_cols<7; ++n_cols)
    {
      deallog << "n_cols = " << n_cols << std::endl;
      ChunkSparsityPattern sp (5,n_cols,3,chunk_size);
      for (unsigned int i=0; i<5; ++i)
        for (unsigned int j=0; j<n_cols; ++j)
          if ((i+2*j+1) % 3 == 0)
            sp.add (i,j);
      sp.compress ();

      ChunkSparseMatrix<double> m(sp);

      // first set a few entries
      for (unsigned int i=0; i<m.m(); ++i)
        for (unsigned int j=0; j<m.n(); ++j)
          if ((i+2*j+1) % 3 == 0)
            m.set (i,j, i*j*.5+.5);

      // next perform a matrix-vector product using the entries as given by
      // the iterator and compare it with the exact value
      Vector<double> src(m.n()), dst(m.m()), dst_ref(m.m());
      for (unsigned int i=0; i<src.size(); ++i)
        src(i) = (double)Testing::rand()/RAND_MAX;
      for (unsigned int i=0; i<m.m(); ++i)
        {
          double sum = 0;
          for (ChunkSparseMatrix<double>::const_iterator it = m.begin(i);
               it != m.end(i); ++it)
            sum += it->value() * src(it->column());
          dst(i) = sum;
        }
      m.vmult(dst_ref, src);
      dst -= dst_ref;
      deallog << "Error in matrix-vector product done via iterator: "
              << dst.linfty_norm() << std::endl;
    }
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
