// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Test that ComputeIndexOwner can handle IndexSets with gaps (also
// in the context of Vectors and Partitioners)

#include <deal.II/base/mpi.h>
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
    const auto [owning_ranks_of_ghosts, import_data] =
      Utilities::MPI::compute_index_owner_and_requesters(local_owned,
                                                         local_relevant,
                                                         comm);

    deallog << "owning_ranks_of_ghosts:" << std::endl;
    for (auto i : owning_ranks_of_ghosts)
      deallog << i << ' ';
    deallog << std::endl;

    deallog << "requesters:" << std::endl;
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
      for (unsigned int j = 0; j < v.ghost_targets()[i].second; ++j)
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
