// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test that it is possible to instantiate a LinearOperator object for all
// different kinds of Trilinos matrices and vectors
// TODO: A bit more tests...

#include <deal.II/lac/linear_operator.h>
#include <deal.II/lac/trilinos_block_sparse_matrix.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_tpetra_sparse_matrix.h>
#include <deal.II/lac/trilinos_tpetra_vector.h>
#include <deal.II/lac/trilinos_vector.h>

#include "../tests.h"


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  initlog();
  deallog << std::setprecision(10);

  {
    TrilinosWrappers::SparseMatrix a;
    auto op_a  = linear_operator<TrilinosWrappers::MPI::Vector>(a);
    auto op_a2 = linear_operator<TrilinosWrappers::MPI::Vector>(a);
  }

#ifdef DEAL_II_TRILINOS_WITH_TPETRA
  {
    LinearAlgebra::TpetraWrappers::SparseMatrix<double> a;
    auto                                                op_a =
      linear_operator<LinearAlgebra::TpetraWrappers::Vector<double>>(a);
    auto op_a2 =
      linear_operator<LinearAlgebra::TpetraWrappers::Vector<double>>(a);
  }
#endif

  {
    TrilinosWrappers::BlockSparseMatrix b;
    auto op_b  = linear_operator<TrilinosWrappers::MPI::BlockVector>(b);
    auto op_b3 = linear_operator<TrilinosWrappers::MPI::BlockVector>(b);
  }

  deallog << "OK" << std::endl;

  return 0;
}
