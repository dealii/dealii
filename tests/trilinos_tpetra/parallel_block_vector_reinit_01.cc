// ---------------------------------------------------------------------
//
// Copyright (C) 2022 - 2023 by the deal.II authors
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


// test various constructors and reinit functions

#include <deal.II/base/index_set.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/trilinos_tpetra_parallel_block_vector.h>

#include <vector>

#include "../tests.h"


void
test()
{
  MPI_Comm comm = MPI_COMM_WORLD;

  const unsigned int myid    = Utilities::MPI::this_mpi_process(comm);
  const unsigned int n_procs = Utilities::MPI::n_mpi_processes(comm);

  constexpr unsigned int n_blocks                     = 3;
  constexpr unsigned int n_indices_per_proc_and_block = 2;

  const unsigned int n_indices_per_block =
    n_indices_per_proc_and_block * n_procs;
  const unsigned int n_indices = n_blocks * n_indices_per_block;

  // set up partitioning
  std::vector<IndexSet> owned_indexsets;
  std::vector<IndexSet> relevant_indexsets;
  std::vector<std::shared_ptr<const Utilities::MPI::Partitioner>> partitioners;
  for (unsigned int b = 0; b < n_blocks; ++b)
    {
      const unsigned int begin = b * n_indices_per_block;

      IndexSet owned(n_indices);
      owned.add_range(begin + n_indices_per_proc_and_block * myid,
                      begin + n_indices_per_proc_and_block * (myid + 1));

      IndexSet relevant(n_indices);
      relevant.add_range(begin + 1, begin + 2);

      partitioners.push_back(
        std::make_shared<const Utilities::MPI::Partitioner>(owned,
                                                            relevant,
                                                            comm));
      owned_indexsets.push_back(std::move(owned));
      relevant_indexsets.push_back(std::move(relevant));
    }

  // create block vectors using different constructors
  LinearAlgebra::TpetraWrappers::MPI::BlockVector block_vector;

  block_vector.reinit(owned_indexsets, comm);
  AssertThrow(block_vector.has_ghost_elements() == false, ExcInternalError());
  deallog << "w/o ghost indices: OK" << std::endl;

  block_vector.reinit(owned_indexsets, relevant_indexsets, comm);
  AssertThrow(block_vector.has_ghost_elements() == true, ExcInternalError());
  deallog << "w/  ghost indices: OK" << std::endl;

  block_vector.reinit(partitioners, /*make_ghosted=*/false);
  AssertThrow(block_vector.has_ghost_elements() == false, ExcInternalError());
  deallog << "partitioners w/o ghost indices: OK" << std::endl;

  block_vector.reinit(partitioners, /*make_ghosted=*/true);
  AssertThrow(block_vector.has_ghost_elements() == true, ExcInternalError());
  deallog << "partitioners w/  ghost indices: OK" << std::endl;
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  const unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  if (myid == 0)
    initlog();

  test();
}
