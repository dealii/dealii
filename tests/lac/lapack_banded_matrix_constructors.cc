// ---------------------------------------------------------------------
//
// Copyright (C) 2018 by the deal.II authors
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

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/lapack_banded_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparse_matrix_ez.h>
#include <deal.II/lac/sparsity_pattern.h>

#include "../tests.h"

// Test the various ways in which one may construct a LAPACKBandedMatrix.

template <typename Number>
void
print_exists(const LAPACKBandedMatrix<Number> &matrix)
{
  for (unsigned int row_n = 0; row_n < matrix.m(); ++row_n)
    {
      for (unsigned int col_n = 0; col_n < matrix.n(); ++col_n)
        {
          deallog << matrix.exists(row_n, col_n) << ' ';
        }
      deallog << std::endl;
    }
}



template <typename Number>
void
print_with_el(const LAPACKBandedMatrix<Number> &matrix)
{
  for (unsigned int row_n = 0; row_n < matrix.n(); ++row_n)
    {
      for (unsigned int col_n = 0; col_n < matrix.n(); ++col_n)
        {
          deallog << std::setw(4) << int(matrix.el(row_n, col_n)) << ' ';
        }
      deallog << std::endl;
    }
}



template <typename Number>
void
print_with_catch(const LAPACKBandedMatrix<Number> &matrix)
{
  // This is slow since we catch a lot of exceptions
#ifdef DEBUG
  deal_II_exceptions::disable_abort_on_exception();
  for (unsigned int row_n = 0; row_n < matrix.m(); ++row_n)
    {
      for (unsigned int col_n = 0; col_n < matrix.n(); ++col_n)
        {
          try
            {
              const int value = int(matrix(row_n, col_n));
              deallog << value << ' ';
            }
          catch (...)
            {
              deallog << "  ";
            }
        }
      deallog << std::endl;
    }
#endif
}

int
main()
{
  initlog();

  {
    // just upper triangular: note that we add superdiagonals as long as at
    // least one row has an entry in the superdiagonal
    const unsigned int n_rows = 10;
    SparsityPattern    sparsity_pattern(n_rows, n_rows);
    for (unsigned int row_n = 0; row_n < n_rows; ++row_n)
      {
        sparsity_pattern.add(row_n, row_n);
        if (row_n < 4)
          sparsity_pattern.add(row_n, row_n + 1);
      }
    sparsity_pattern.compress();

    LAPACKBandedMatrix<float> banded(sparsity_pattern);
    print_exists(banded);
  }


  deallog << std::endl << "SparsityPattern" << std::endl;
  {
    // three lower diagonals, one upper diagonal
    const unsigned int n_rows = 10;
    SparsityPattern    sparsity_pattern(n_rows, n_rows);
    for (unsigned int row_n = 0; row_n < n_rows; ++row_n)
      {
        sparsity_pattern.add(row_n, row_n);
        if (row_n < n_rows - 1)
          sparsity_pattern.add(row_n, row_n + 1);
        if (0 < row_n)
          sparsity_pattern.add(row_n, row_n - 1);
      }
    sparsity_pattern.compress();

    LAPACKBandedMatrix<float> banded(sparsity_pattern);
    print_exists(banded);
  }

  deallog << std::endl << "SparseMatrix<float>" << std::endl;
  {
    // three lower diagonals, one upper diagonal
    const unsigned int n_rows = 10;
    SparsityPattern    sparsity_pattern(n_rows, n_rows);
    for (unsigned int row_n = 0; row_n < n_rows; ++row_n)
      {
        sparsity_pattern.add(row_n, row_n);
        if (row_n < n_rows - 1)
          sparsity_pattern.add(row_n, row_n + 1);
        if (0 < row_n)
          sparsity_pattern.add(row_n, row_n - 1);
      }
    sparsity_pattern.compress();

    SparseMatrix<float> mat(sparsity_pattern);
    for (unsigned int row_n = 0; row_n < n_rows; ++row_n)
      {
        mat.set(row_n, row_n, 2.0);
        if (0 < row_n)
          mat.set(row_n, row_n - 1, 1.0);
        if (row_n < n_rows - 1)
          mat.set(row_n, row_n + 1, 1.0);
      }

    LAPACKBandedMatrix<float> banded(mat);
    print_with_el(banded);
    deallog << std::endl;
    print_with_catch(banded);
  }

  deallog << std::endl << "SparseMatrixEZ<double>" << std::endl;
  {
    const int              n_rows = 5;
    SparseMatrixEZ<double> mat(n_rows, n_rows);
    for (int row_n = 0; row_n < n_rows; ++row_n)
      {
        mat.set(row_n, row_n, 2.0 + row_n);
        if (0 < row_n)
          mat.set(row_n, row_n - 1, 1.0 - row_n);
        if (row_n < n_rows - 1)
          mat.set(row_n, row_n + 1, 1.0);
      }

    LAPACKBandedMatrix<double> banded(mat);
    print_with_el(banded);
    deallog << std::endl;
  }

  deallog << std::endl << "FullMatrix<double>" << std::endl;
  {
    const int          n_rows = 5;
    FullMatrix<double> mat(n_rows, n_rows);
    for (int row_n = 0; row_n < n_rows; ++row_n)
      {
        mat.set(row_n, row_n, 2.0 + row_n);
        if (0 < row_n)
          mat.set(row_n, row_n - 1, 1.0);
        if (row_n < n_rows - 1)
          mat.set(row_n, row_n + 1, 1.0);
      }

    // verify that we did not add extra diagonals
    LAPACKBandedMatrix<double> banded(mat);
    print_with_el(banded);
    deallog << std::endl;
    print_with_catch(banded);
  }
}
