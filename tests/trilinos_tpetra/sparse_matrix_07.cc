// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2024 - 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------



// like sparse_matrix_06, but in parallel

#include <deal.II/base/utilities.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/trilinos_tpetra_sparse_matrix.h>

#include <iostream>

#include "../tests.h"


int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());

  MPILogInitAll log;

  int      my_rank = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  IndexSet locally_owned_dofs(5);
  if (my_rank == 0)
    locally_owned_dofs.add_range(0, 3);
  else
    locally_owned_dofs.add_range(3, 5);

  IndexSet locally_relevant_dofs = locally_owned_dofs;

  DynamicSparsityPattern dsp(locally_relevant_dofs);
  if (my_rank == 0)
    {
      dsp.add(1, 2);
      dsp.add(2, 3);
    }
  else
    {
      dsp.add(3, 4);
      dsp.add(4, 3);
    }

  SparsityTools::distribute_sparsity_pattern(dsp,
                                             locally_owned_dofs,
                                             MPI_COMM_WORLD,
                                             locally_relevant_dofs);

  LinearAlgebra::TpetraWrappers::SparseMatrix<double, MemorySpace::Default> matrix;
  matrix.reinit(locally_owned_dofs, dsp, MPI_COMM_WORLD);

  deallog << "Original:" << std::endl;
  matrix.print(deallog.get_file_stream());

  // now reinitialize a second Trilinos matrix
  LinearAlgebra::TpetraWrappers::SparseMatrix<double, MemorySpace::Default>
    tmatrix;
  tmatrix.reinit(matrix);

  deallog << "Copy:" << std::endl;
  tmatrix.print(deallog.get_file_stream());
  deallog << "OK" << std::endl;
}
