// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// test LinearAlgebra::TpetraWrappers::SparseMatrix<double>::reinit and add
// when one MPI process doesn't own any rows.

#include <deal.II/base/utilities.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/trilinos_tpetra_sparse_matrix.h>

#include <iostream>

#include "../tests.h"


int
main(int argc, char **argv)
{
  initlog();

  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());

  int      my_rank = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  IndexSet locally_owned_dofs(10);
  if (my_rank == 0)
    locally_owned_dofs.add_range(0, 10);
  IndexSet locally_relevant_dofs = locally_owned_dofs;

  DynamicSparsityPattern dsp(locally_relevant_dofs);
  if (my_rank == 0)
    dsp.add(0, 0);
  SparsityTools::distribute_sparsity_pattern(dsp,
                                             locally_owned_dofs,
                                             MPI_COMM_WORLD,
                                             locally_relevant_dofs);

  LinearAlgebra::TpetraWrappers::SparseMatrix<double> matrix;
  matrix.reinit(locally_owned_dofs, dsp, MPI_COMM_WORLD);

  Assert(matrix.is_compressed() == true, ExcInternalError());

  if (my_rank == 0)
    {
      matrix.set(0, 0, 10);
      Assert(matrix.is_compressed() == false, ExcInternalError());
    }
  else
    Assert(matrix.is_compressed() == true, ExcInternalError());

  matrix.compress(VectorOperation::insert);
  deallog << "OK" << std::endl;
}
