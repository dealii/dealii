// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/base/exceptions.h>

#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparse_matrix_ez.h>
#include <deal.II/lac/vector.h>

#include <iostream>

#include "../tests.h"


/**
 * Finite difference diagonal matrix on uniform grid.
 */

class FDDiagMatrix
{
public:
  /**
   * Constructor specifying grid resolution.
   */
  FDDiagMatrix(unsigned int nx, unsigned int ny);

  /**
   * Generate the matrix structure.
   */
  template <typename SP>
  void
  diag_structure(SP &structure) const;

  /**
   * Fill the matrix with values.
   */
  template <typename MatrixType>
  void
  diag(MatrixType &) const;

  template <typename number>
  void
  gnuplot_print(std::ostream &, const Vector<number> &) const;

private:
  /**
   * Number of gridpoints in x-direction.
   */
  unsigned int nx;

  /**
   * Number of gridpoints in y-direction.
   */
  unsigned int ny;
};


// --------------- inline and template functions -----------------

inline FDDiagMatrix::FDDiagMatrix(unsigned int nx, unsigned int ny)
  : nx(nx)
  , ny(ny)
{}



template <typename SP>
inline void
FDDiagMatrix::diag_structure(SP &structure) const
{
  for (unsigned int i = 0; i <= ny - 2; ++i)
    {
      for (unsigned int j = 0; j <= nx - 2; ++j)
        {
          // Number of the row to be entered
          unsigned int row = j + (nx - 1) * i;
          structure.add(row, row);
        }
    }
}


template <typename MatrixType>
inline void
FDDiagMatrix::diag(MatrixType &A) const
{
  for (unsigned int i = 0; i <= ny - 2; ++i)
    {
      for (unsigned int j = 0; j <= nx - 2; ++j)
        {
          // Number of the row to be entered
          unsigned int row = j + (nx - 1) * i;
          A.set(row, row, 2. / std::log(2.0 + row));
        }
    }
}

template <typename number>
inline void
FDDiagMatrix::gnuplot_print(std::ostream &s, const Vector<number> &V) const
{
  for (unsigned int i = 0; i <= ny - 2; ++i)
    {
      for (unsigned int j = 0; j <= nx - 2; ++j)
        {
          // Number of the row to be entered
          unsigned int row = j + (nx - 1) * i;
          s << (j + 1) / (float)nx << '\t' << (i + 1) / (float)ny << '\t'
            << V(row) << std::endl;
        }
      s << std::endl;
    }
  s << std::endl;
}


/**
 * Finite difference matrix for laplace operator in 1D on uniform grid.
 */

class FD1DLaplaceMatrix
{
public:
  /**
   * Constructor specifying grid resolution.
   */
  FD1DLaplaceMatrix(unsigned int n);

  /**
   * Generate the matrix structure.
   */
  template <typename SP>
  void
  three_point_structure(SP &structure) const;

  /**
   * Fill the matrix with values.
   */
  template <typename MatrixType>
  void
  three_point(MatrixType &) const;

private:
  /**
   * Number of gridpoints in x-direction.
   */
  unsigned int n;
};


// --------------- inline and template functions -----------------

inline FD1DLaplaceMatrix::FD1DLaplaceMatrix(unsigned int n)
  : n(n)
{}



template <typename SP>
inline void
FD1DLaplaceMatrix::three_point_structure(SP &structure) const
{
  for (unsigned int i = 0; i <= n - 2; ++i)
    {
      structure.add(i, i);
      if (i > 0)
        structure.add(i, i - 1);
      if (i < n - 1)
        structure.add(i, i + 1);
    }
}


template <typename MatrixType>
inline void
FD1DLaplaceMatrix::three_point(MatrixType &A) const
{
  for (unsigned int i = 0; i <= n - 2; ++i)
    {
      A.set(i, i, 2.0);
      if (i > 0)
        A.set(i, i - 1, -1.0);

      if (i < n - 2)
        A.set(i, i + 1, -1.0);
    }
}
