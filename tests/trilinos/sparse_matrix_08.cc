// ---------------------------------------------------------------------
//
// Copyright (C) 2017 by the deal.II authors
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


// Test move constructor of sparse matrix and check that the matrix that has
// been move can be reinitialized

#include <deal.II/base/utilities.h>

#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_sparsity_pattern.h>

#include <iostream>

#include "../tests.h"

int
main(int argc, char **argv)
{
  initlog();

  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());

  IndexSet partitioning(3);

  partitioning.add_range(0, 3);

  // Add element (2,1) to the matrix
  TrilinosWrappers::SparsityPattern sp(partitioning);
  sp.add(2, 1);
  sp.compress();

  TrilinosWrappers::SparseMatrix A(sp);
  A.add(2, 1, 2.0);
  A.compress(VectorOperation::add);

  // Check that the B stole the data of A
  TrilinosWrappers::SparseMatrix B(std::move(A));
  for (unsigned int i = 0; i < 3; ++i)
    for (unsigned int j = 0; j < 3; ++j)
      {
        if ((i == 2) && (j == 1))
          {
            AssertThrow(B.el(i, j) == 2, ExcInternalError());
          }
        else
          {
            AssertThrow(B.el(i, j) == 0, ExcInternalError());
          }
      }

  // Check that A can be reinitialized
  TrilinosWrappers::SparsityPattern sp_2(partitioning);
  sp_2.add(1, 2);
  sp_2.compress();
  A.reinit(sp_2);
  A.add(1, 2, 2.0);
  A.compress(VectorOperation::add);
  for (unsigned int i = 0; i < 3; ++i)
    for (unsigned int j = 0; j < 3; ++j)
      {
        if ((i == 1) && (j == 2))
          {
            AssertThrow(A.el(i, j) == 2, ExcInternalError());
          }
        else
          {
            AssertThrow(A.el(i, j) == 0, ExcInternalError());
          }
      }

  deallog << "OK" << std::endl;
}
