// ---------------------------------------------------------------------
//
// Copyright (C) 2015 by the deal.II authors
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

// Tests for the LinearOperator template with
//   dealii::BlockVector<double>
//   dealii::BlockSparseMatrix<double>

#include "../tests.h"

#include <deal.II/lac/linear_operator.h>
#include <deal.II/lac/trilinos_block_sparse_matrix.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_vector.h>

using namespace dealii;

int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  initlog();
  deallog << std::setprecision(10);

  TrilinosWrappers::SparseMatrix a;
  auto op_a  = linop<TrilinosWrappers::MPI::Vector>(a);

  TrilinosWrappers::BlockSparseMatrix b;
  auto op_b = linop<TrilinosWrappers::MPI::BlockVector>(b);

  static const int dim = 2;

  return 0;
}
