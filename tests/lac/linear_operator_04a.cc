// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2018 by the deal.II authors
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

// Test that it is possible to use a PackagedOperation created by
// operator*() of a LinearOperator object for Trilinos matrices and vectors

#include <deal.II/lac/linear_operator_tools.h>
#include <deal.II/lac/trilinos_block_sparse_matrix.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_vector.h>

#include "../tests.h"


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  initlog();
  deallog << std::setprecision(10);

  typedef TrilinosWrappers::MPI::Vector  vector_t;
  typedef TrilinosWrappers::SparseMatrix matrix_t;

  matrix_t a(5U, 5U, 3U);
  a.compress(VectorOperation::add);

  auto     op_a = linear_operator<vector_t>(a);
  vector_t u, res;
  op_a.reinit_domain_vector(u, false);
  res = op_a * u;
  // ^^ this was not working, whereas op_a.vmult(res,u) did.

  deallog << "OK" << std::endl;

  return 0;
}
