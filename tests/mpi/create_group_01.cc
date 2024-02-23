// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check that Utilities::MPI::create_group is only collective over the group of
// processes that actually want to create the communicator.
// In p4est-2.0, sc_init() leads to a deadlock when MPI_Comm_create_group is
// called (see https://github.com/cburstedde/p4est/issues/30). Hence, this test
// also ensures that the treatment of sc_init() is appropriate.

#include "../tests.h"

int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv, 1);
  MPILogInitAll                    log_all;

  const MPI_Comm global_comm = MPI_COMM_WORLD;

  const unsigned int my_id   = Utilities::MPI::this_mpi_process(global_comm);
  const unsigned int n_ranks = Utilities::MPI::n_mpi_processes(global_comm);

  // Get the group of processes in MPI_COMM_WORLD
  MPI_Group world_group;
  int       ierr = MPI_Comm_group(MPI_COMM_WORLD, &world_group);
  AssertThrowMPI(ierr);

  const int        n = n_ranks / 2;
  std::vector<int> ranks;
  for (unsigned int i = 0; i < n_ranks; i += 2)
    ranks.push_back(i);

  // Construct a group containing all of the even ranks in world_group
  MPI_Group even_group;
  ierr = MPI_Group_incl(world_group, n, ranks.data(), &even_group);
  AssertThrowMPI(ierr);

  if (my_id % 2 == 0)
    {
      MPI_Comm group_comm;
      ierr =
        Utilities::MPI::create_group(global_comm, even_group, 0, &group_comm);
      AssertThrowMPI(ierr);
    }
  deallog << "OK" << std::endl;
}
