// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// test MPI::Partitioner setup for the case the owned indices are not onto,
// i.e., we skip some index (related test: parallel_vector_16 and
// compute_index_owner_01)

#include <deal.II/base/partitioner.h>

#include "../tests.h"

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize init(argc, argv, 1);

  MPILogInitAll log;

  IndexSet owned(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD) + 1);
  owned.add_index(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) + 1);

  IndexSet ghosted(owned.size());
  ghosted.add_index(1 + (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) + 1) %
                          Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD));

  Utilities::MPI::Partitioner part(owned, ghosted, MPI_COMM_WORLD);

  deallog << "ghost targets: ";
  for (const auto &p : part.ghost_targets())
    deallog << 'p' << p.first << " n_indices=" << p.second << "   ";
  deallog << std::endl;
  deallog << "import targets: ";
  for (const auto &p : part.import_targets())
    deallog << 'p' << p.first << " n_indices=" << p.second << "   ";
  deallog << std::endl;
  deallog << "import indices: ";
  for (const auto &p : part.import_indices())
    deallog << '[' << p.first << ' ' << p.second << ") ";
  deallog << std::endl;
}
