// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2011 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// TrilinosWrappers::SparseMatrix::el() and operator() had problems
// looking up indices in rectangular matrices because they
// accidentally used the row instead of the column map. verify that
// this is now fixed
//
// this testcase is reduced from one contributed by Habib Talavatifard

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

  IndexSet row_partitioning(3);
  IndexSet col_partitioning(4);

  row_partitioning.add_range(0, 3);
  col_partitioning.add_range(0, 4);

  // Add element (2,3) to the matrix
  TrilinosWrappers::SparsityPattern sp(row_partitioning, col_partitioning);
  sp.add(2, 3);
  sp.compress();

  TrilinosWrappers::SparseMatrix A(sp);
  A.add(2, 3, 2.0);
  A.compress(VectorOperation::add);

  // verify that entry (2,3) is
  // indeed what we expect. verify
  // that both methods of accessing
  // the entry work
  AssertThrow(A.el(2, 3) == 2, ExcInternalError());
  AssertThrow(A(2, 3) == 2, ExcInternalError());

  deallog << "OK" << std::endl;
}
