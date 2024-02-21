// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2004 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// compare collective adding of elements in a trilinos matrix using
// TrilinosWrappers::SparseMatrix::add() with element-wise setting

#include <deal.II/base/utilities.h>

#include <deal.II/lac/trilinos_sparse_matrix.h>

#include <iostream>

#include "../tests.h"


void
test(TrilinosWrappers::SparseMatrix &m)
{
  TrilinosWrappers::SparseMatrix m2(m.m(), m.n(), m.m() / 3 + 1);

  // first set a few entries one-by-one
  for (unsigned int i = 0; i < m.m(); ++i)
    for (unsigned int j = 0; j < m.n(); ++j)
      if ((i + 2 * j + 1) % 3 == 0)
        {
          m.set(i, j, i * j * .5 + .5);
          m2.set(i, j, 0.);
        }

  m.compress(VectorOperation::insert);
  m2.compress(VectorOperation::insert);

  // now add the same elements row-wise
  {
    std::vector<types::global_dof_index> col_indices(m.n() / 3 + 1);
    std::vector<double>                  col_values(m.n() / 3 + 1);
    for (unsigned int i = 0; i < m.m(); ++i)
      {
        unsigned int col_index = 0;
        // count the number of elements in this
        // row
        for (unsigned int j = 0; j < m.n(); ++j)
          if ((i + 2 * j + 1) % 3 == 0)
            ++col_index;

        col_indices.resize(col_index);
        col_values.resize(col_index);
        col_index = 0;

        // extract column values
        for (unsigned int j = 0; j < m.n(); ++j)
          if ((i + 2 * j + 1) % 3 == 0)
            {
              col_indices[col_index] = j;
              col_values[col_index]  = i * j * .5 + .5;
              col_index++;
            }

        m2.add(i, col_indices, col_values);
      }
  }

  m2.compress(VectorOperation::add);

  // subtract the matrix m from this one,
  // we should get a zero matrix
  m2.add(-1.0, m);
  // calculate the Frobenius norm of the
  // matrix in order to check whether all
  // elements really are zero
  double norm = m2.frobenius_norm();
  AssertThrow(norm == 0, ExcInternalError());

  deallog << "OK" << std::endl;
}



int
main(int argc, char **argv)
{
  initlog();

  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());

  try
    {
      {
        TrilinosWrappers::SparseMatrix m(5U, 5U, 3U);

        test(m);
      }
    }
  catch (const std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
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
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    };
}
