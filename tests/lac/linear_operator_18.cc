// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test that it is possible to use a PackagedOperation created by
// operator*() of a LinearOperator object for Trilinos matrices and vectors

#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/linear_operator_tools.h>
#include <deal.II/lac/trilinos_block_sparse_matrix.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>

#include "../tests.h"


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  initlog();
  deallog << std::setprecision(10);

  using vector_t = LinearAlgebra::distributed::Vector<double>;
  using matrix_t = TrilinosWrappers::SparseMatrix;

  matrix_t a(5U, 4U, 3U);
  a.compress(VectorOperation::add);

  auto     op_a = linear_operator<vector_t>(a);
  vector_t u, res;
  op_a.reinit_domain_vector(u, false);
  res = op_a * u;

  auto op_a_transpose = transpose_operator(op_a);
  op_a.reinit_range_vector(u, false);
  res = op_a_transpose * u;

  deallog << "OK" << std::endl;

  return 0;
}
