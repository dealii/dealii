// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2020 by the deal.II authors
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



// Test that ComputeIndexOwner can handle IndexSets with gaps (also
// in the context of Vectors and Partitioners)

#include <deal.II/base/mpi.h>
#include <deal.II/base/mpi_compute_index_owner_internal.h>
#include <deal.II/base/partitioner.h>

#include "../tests.h"


void
test()
{
  MPI_Comm comm = MPI_COMM_WORLD;


  unsigned int myid    = Utilities::MPI::this_mpi_process(comm);
  unsigned int numproc = Utilities::MPI::n_mpi_processes(comm);

  unsigned int offset = 1;
  unsigned int size   = 1;

  const unsigned int global_size = (numproc + offset) * size;

  // locals
  IndexSet local_owned(global_size);
  local_owned.add_range((myid + offset) * size, (myid + offset + 1) * size);

  // ghosts
  IndexSet local_relevant(global_size);
  local_relevant.add_range(((myid + 1) % numproc + offset) * size,
                           ((myid + 1) % numproc + offset + 1) * size);

  deallog << "local" << std::endl;
  local_owned.print(deallog);
  deallog << "ghost" << std::endl;
  local_relevant.print(deallog);

  {
    std::vector<unsigned int> owning_ranks_of_ghosts(
      local_relevant.n_elements());

    Utilities::MPI::internal::ComputeIndexOwner::ConsensusAlgorithmsPayload
      process(local_owned, local_relevant, comm, owning_ranks_of_ghosts, true);

    Utilities::MPI::ConsensusAlgorithms::Selector<
      std::pair<types::global_dof_index, types::global_dof_index>,
      unsigned int>
      consensus_algorithm(process, comm);
    consensus_algorithm.run();

    deallog << "owning_ranks_of_ghosts:" << std::endl;
    for (auto i : owning_ranks_of_ghosts)
      deallog << i << " ";
    deallog << std::endl;

    deallog << "requesters:" << std::endl;
    std::map<unsigned int, IndexSet> import_data = process.get_requesters();
    for (const auto &m : import_data)
      {
        deallog << m.first << ": ";
        m.second.print(deallog);
        deallog << std::endl;
      }
  }


  {
    Utilities::MPI::Partitioner v(local_owned, local_relevant, comm);

    for (unsigned int i = 0; i < v.ghost_targets().size(); ++i)
      for (unsigned int j = 0; j < v.ghost_targets()[i].second; j++)
        deallog << v.ghost_targets()[i].first << std::endl;
  }
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    all;

  test();
}
