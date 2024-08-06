// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------



// compare collective setting of elements in a trilinos matrix using
// LinearAlgebra::TpetraWrappers::SparseMatrix<double>::set() with element-wise
// setting

#include <deal.II/base/utilities.h>

#include <deal.II/lac/trilinos_tpetra_sparse_matrix.h>

#include <iostream>

#include "../tests.h"


void
test(
  LinearAlgebra::TpetraWrappers::SparseMatrix<double, MemorySpace::Default> &m)
{
  // first set a few entries one-by-one
  for (unsigned int i = 0; i < m.m(); ++i)
    for (unsigned int j = 0; j < m.n(); ++j)
      if ((i + 2 * j + 1) % 3 == 0)
        m.set(i, j, i * j * .5 + .5);

  m.compress(VectorOperation::insert);

  LinearAlgebra::TpetraWrappers::SparseMatrix<double, MemorySpace::Default> m2(
    m.m(), m.n(), m.n() / 3 + 1);

  // now set the same elements row-wise
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

        m2.set(i, col_indices, col_values);
      }
  }

  m2.compress(VectorOperation::insert);

  // subtract the matrix m from this one,
  // we should get a zero matrix.
  m2.add(-1.0, m);
  // Calculate the Frobenius norm of the
  // matrix in order to check whether all
  // elements really are zero.
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
        LinearAlgebra::TpetraWrappers::SparseMatrix<double,
                                                    MemorySpace::Default>
          m(5U, 5U, 3U);
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
