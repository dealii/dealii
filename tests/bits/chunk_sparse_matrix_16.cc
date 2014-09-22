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



// compare matrix-vector product in chunk sparse matrix with usual vmult and
// extract_row_copy.

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
      // the entries in extract_row_copy and compare it with the exact value
      Vector<double> src(m.n()), dst(m.m()), dst_ref(m.m());
      for (unsigned int i=0; i<src.size(); ++i)
        src(i) = (double)Testing::rand()/RAND_MAX;
      std::vector<types::global_dof_index> indices(sp.max_entries_per_row());
      std::vector<double> values(sp.max_entries_per_row());
      for (unsigned int i=0; i<m.m(); ++i)
        {
          types::global_dof_index n_entries = numbers::invalid_unsigned_int;
          m.extract_row_copy(i, values.size(), n_entries, &indices[0], &values[0]);
          double sum = 0;
          for (unsigned int j=0; j<n_entries; ++j)
            sum += values[j] * src(indices[j]);
          dst(i) = sum;
        }
      m.vmult(dst_ref, src);
      dst -= dst_ref;
      deallog << "Error in matrix-vector product done via extract_row_copy: "
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
