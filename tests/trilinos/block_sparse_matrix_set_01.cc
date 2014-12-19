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



// compare collective setting elements in a trilinos matrix using
// TrilinosWrappers::BlockSparseMatrix::set() and a FullMatrix<double> with
// setting the same elements on an entry-by-entry. Use the entries as they
// would result from a vector-valued problem with a Laplace operator in 1D
// in each component

#include "../tests.h"
#include <deal.II/base/utilities.h>
#include <deal.II/lac/trilinos_block_sparse_matrix.h>
#include <deal.II/lac/full_matrix.h>
#include <fstream>
#include <iostream>


void test ()
{
  // first set a few entries one-by-one in
  // a small matrix
  const unsigned int block_size = 16;
  TrilinosWrappers::SparseMatrix m_small (block_size,block_size,3U);
  for (unsigned int i=0; i<block_size; ++i)
    for (unsigned int j=0; j<block_size; ++j)
      if (std::fabs((double)i-j) < 2)
        {
          double value;
          if (i == j)
            value = 1.;      // for this example, set the values to
          // one and not some to two, since we do
          // not insert elements into the sparsity
          // pattern, but only set already present
          // ones.
          else
            value = -1.;

          m_small.set (i,j, value);
        }
  m_small.compress();

  // Then build two matrices consisting of
  // several copies of the small
  // matrix. The matrix m2 will be filled
  // later.
  TrilinosWrappers::BlockSparseMatrix m, m2;
  m.reinit(3,2);
  m2.reinit(3,2);
  for (unsigned int block_row = 0; block_row<m.n_block_rows(); ++block_row)
    for (unsigned int block_col = 0; block_col<m.n_block_cols(); ++block_col)
      m.block(block_row, block_col).copy_from(m_small);
  m.collect_sizes();

  // fill the second matrix with the
  // sparsity pattern (but not the
  // values)
  for (unsigned int block_row = 0; block_row<m.n_block_rows(); ++block_row)
    for (unsigned int block_col = 0; block_col<m.n_block_cols(); ++block_col)
      m2.block(block_row, block_col).reinit(m_small);
  m2.collect_sizes();

  // now add the same elements from a full
  // matrix
  {
    FullMatrix<double> full_matrix(2*m.n_block_rows(),2*m.n_block_cols());
    for (unsigned int block_row = 0; block_row<m.n_block_rows(); ++block_row)
      for (unsigned int block_col = 0; block_col<m.n_block_cols(); ++block_col)
        {
          full_matrix(block_row*2,block_col*2) = 1;
          full_matrix(block_row*2+1,block_col*2+1) = 1.;
          full_matrix(block_row*2,block_col*2+1) = -1;
          full_matrix(block_row*2+1,block_col*2) = -1.;
        }

    std::vector<types::global_dof_index> local_row_indices (2*m.n_block_rows());
    std::vector<types::global_dof_index> local_col_indices (2*m.n_block_cols());

    for (unsigned int i=0; i<block_size-1; ++i)
      {
        for (unsigned int block_row = 0; block_row<m.n_block_rows(); ++block_row)
          {
            local_row_indices[2*block_row] = block_row*block_size + i;
            local_row_indices[2*block_row+1] = block_row*block_size + i + 1;
          }
        for (unsigned int block_col = 0; block_col<m.n_block_cols(); ++block_col)
          {
            local_col_indices[2*block_col] = block_col*block_size + i;
            local_col_indices[2*block_col+1] = block_col*block_size + i + 1;
          }

        m2.set (local_row_indices, local_col_indices, full_matrix);
      }
  }

  m2.compress();

  // subtract the matrix m from this one,
  // we should get a zero matrix
  double norm = 0;
  for (unsigned int block_row = 0; block_row<m.n_block_rows(); ++block_row)
    for (unsigned int block_col = 0; block_col<m.n_block_cols(); ++block_col)
      {
        m2.block(block_row,block_col).add(-1.0, m.block(block_row,block_col));

        // calculate the Frobenius norm of the
        // matrix in order to check whether all
        // elements really are zero
        norm += m2.block(block_row,block_col).frobenius_norm();
      }
  Assert (norm == 0, ExcInternalError());

  deallog << "OK" << std::endl;
}



int main (int argc,char **argv)
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

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
}
