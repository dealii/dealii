// ---------------------------------------------------------------------
//
// Copyright (C) 2016 by the deal.II authors
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

// Test that it is possible to use a operator*() of LinearOperator object for
// Trilinos matrices and vectors

#include "../tests.h"

#include <deal.II/lac/linear_operator_tools.h>
#include <deal.II/lac/trilinos_block_sparse_matrix.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_vector.h>

using namespace dealii;

int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv, 1);

  initlog();
  deallog << std::setprecision(10);

  typedef TrilinosWrappers::MPI::Vector vector_t;
  typedef TrilinosWrappers::SparseMatrix matrix_t;

  matrix_t a;

  auto op_a  = linear_operator<vector_t>(a);
  vector_t u,res;
  res = op_a * u;
  // ^^ this was not working, whereas
  // op_a.vmult(res,u) did.

  deallog << "OK" << std::endl;

  return 0;
}
