// ---------------------------------------------------------------------
//
// Copyright (C) 2019 by the deal.II authors
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
