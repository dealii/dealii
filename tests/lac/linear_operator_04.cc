// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2018 by the deal.II authors
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

// Test that it is possible to instantiate a LinearOperator object for all
// different kinds of Trilinos matrices and vectors
// TODO: A bit more tests...

#include <deal.II/lac/linear_operator.h>
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

  TrilinosWrappers::SparseMatrix a;

  auto op_a  = linear_operator<TrilinosWrappers::MPI::Vector>(a);
  auto op_a2 = linear_operator<TrilinosWrappers::MPI::Vector>(a);

  TrilinosWrappers::BlockSparseMatrix b;

  auto op_b  = linear_operator<TrilinosWrappers::MPI::BlockVector>(b);
  auto op_b3 = linear_operator<TrilinosWrappers::MPI::BlockVector>(b);

  deallog << "OK" << std::endl;

  return 0;
}
