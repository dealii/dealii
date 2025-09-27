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
// LinearAlgebra::TpetraWrappers::SparseMatrix<double>::set() and a
// FullMatrix<double> with setting the same elements on an element-by-element
// level. Use the entries as they would result from a Laplace operator in 1D.

#include <deal.II/base/utilities.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/trilinos_tpetra_sparse_matrix.h>

#include <iostream>

#include "../tests.h"


void
test(
  LinearAlgebra::TpetraWrappers::SparseMatrix<double, MemorySpace::Default> &m)
{
  LinearAlgebra::TpetraWrappers::SparseMatrix<double, MemorySpace::Default> m2(
    m.m(), m.n(), 3U);

  // first set a few entries one-by-one and
  // initialize the sparsity pattern for m2
  for (unsigned int i = 0; i < m.m(); ++i)
    for (unsigned int j = 0; j < m.n(); ++j)
      if (std::fabs((double)i - j) < 2)
        {
          double value;
          if (i == j)
            if (i > 0 && i < m.m() - 1)
              value = 2.;
            else
              value = 1.;
          else
            value = -1.;

          m.set(i, j, value);
        }
  m.compress(VectorOperation::insert);
  // now add the same elements from a full
  // matrix (if the same element is set
  // more than once, Trilinos adds them)
  {
    FullMatrix<double> full_matrix(2, 2);
    full_matrix(0, 0) = full_matrix(1, 1) = 1.;
    full_matrix(0, 1) = full_matrix(1, 0) = -1.;
    std::vector<types::global_dof_index> local_indices(2);

    for (unsigned int i = 0; i < m.m() - 1; ++i)
      {
        local_indices[0] = i;
        local_indices[1] = i + 1;

        m2.set(local_indices, local_indices, full_matrix);
      }
  }

  m2.compress(VectorOperation::insert);

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
        LinearAlgebra::TpetraWrappers::SparseMatrix<double,
                                                    MemorySpace::Default>
          m(16U, 16U, 3U);
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
