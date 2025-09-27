// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// make sure that block vector copies correctly to a serial vector

#include <deal.II/base/utilities.h>

#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>

#include <iostream>

#include "../tests.h"


int
main(int argc, char **argv)
{
  using BlockVectorLocal = BlockVector<double>;
  using BlockVector      = TrilinosWrappers::MPI::BlockVector;

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;

  const MPI_Comm     mpi_communicator = MPI_COMM_WORLD;
  const unsigned int this_mpi_process =
    Utilities::MPI::this_mpi_process(mpi_communicator);
  const unsigned int n_mpi_processes =
    Utilities::MPI::n_mpi_processes(mpi_communicator);

  const unsigned int n_blocks         = 2;
  const unsigned int n_dofs_per_block = 10;

  std::vector<IndexSet> locally_owned_partitioning(n_blocks);
  for (unsigned int b = 0; b < n_blocks; ++b)
    {
      locally_owned_partitioning[b].set_size(n_dofs_per_block);
      locally_owned_partitioning[b].add_range(this_mpi_process *
                                                (n_dofs_per_block / 2),
                                              (this_mpi_process + 1) *
                                                (n_dofs_per_block / 2));
    }

  BlockVector parallel_vector;
  parallel_vector.reinit(locally_owned_partitioning);
  deallog << "Locally owned indices" << std::endl;
  parallel_vector.locally_owned_elements().print(deallog.get_file_stream());

  // Set entries in parallel vector
  for (auto idx : parallel_vector.locally_owned_elements())
    {
      parallel_vector[idx] = 10.0 * idx;
    }
  deallog << "Parallel vector" << std::endl;
  parallel_vector.print(deallog.get_file_stream());

  // Copy distributed vector to local vector
  deallog << "Localized vector (copy constructor)" << std::endl;
  const BlockVectorLocal local_vector_1(parallel_vector);
  local_vector_1.print(deallog.get_file_stream());

  // Copy distributed vector to local vector
  deallog << "Localized vector (operator =)" << std::endl;
  const BlockVectorLocal local_vector_2 = parallel_vector;
  local_vector_2.print(deallog.get_file_stream());

  Assert(local_vector_1 == local_vector_2, ExcMessage("Vectors don't match"));

  deallog << "OK" << std::endl;
}
