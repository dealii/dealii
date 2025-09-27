// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check TrilinosWrappers::SparsityPattern::row_length
#include <deal.II/base/index_set.h>
#include <deal.II/base/mpi.h>

#include <deal.II/lac/block_sparsity_pattern.h>

#include "../tests.h"

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  MPILogInitAll mpi_log;

  const unsigned int n = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
  IndexSet           owned_set(n);
  owned_set.add_index(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD));
  IndexSet relevant_set(n);
  relevant_set.add_range(0, n);
  IndexSet                          empty_set;
  TrilinosWrappers::SparsityPattern sparsity(owned_set,
                                             empty_set,
                                             relevant_set,
                                             MPI_COMM_WORLD);

  deallog << "\nBefore compressing sparsity pattern..." << std::endl;
  sparsity.compress();
  deallog << "Compressing sparsity pattern worked!\n\n" << std::endl;

  return 0;
}
