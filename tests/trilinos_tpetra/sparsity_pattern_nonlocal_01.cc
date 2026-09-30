// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2026 - 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------



// Check export of explicitly writable nonlocal rows, including empty patterns.

#include <deal.II/lac/trilinos_tpetra_sparsity_pattern.h>

#include "../tests.h"

void
check(const unsigned int n_columns, const bool insert_entries)
{
  const MPI_Comm     comm    = MPI_COMM_WORLD;
  const unsigned int rank    = Utilities::MPI::this_mpi_process(comm);
  const unsigned int n_ranks = Utilities::MPI::n_mpi_processes(comm);

  // Locally owned rows and columns
  IndexSet locally_owned(n_ranks);
  locally_owned.add_index(rank);
  IndexSet columns(n_columns);
  for (unsigned int c = rank; c < n_columns; c += n_ranks)
    columns.add_index(c);

  // Make all rows writable, even nonlocal rows
  IndexSet writable(n_ranks);
  writable.add_range(0, n_ranks);

  LinearAlgebra::TpetraWrappers::SparsityPattern<MemorySpace::Host> pattern;
  pattern.reinit(locally_owned, columns, writable, comm, n_columns);

  if (insert_entries)
    {
      // Every entry targets a remotely owned row and column.
      const unsigned int entry = (rank + 1) % n_ranks;
      pattern.add(entry, entry);
    }
  pattern.compress();

  AssertThrow(pattern.row_length(rank) == (insert_entries ? 1 : 0),
              ExcInternalError());

  if (insert_entries)
    AssertThrow(pattern.exists(rank, rank), ExcInternalError());

  AssertThrow(pattern.n_nonzero_elements() == (insert_entries ? n_ranks : 0),
              ExcInternalError());
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());

  initlog();

  const unsigned int n_ranks = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);

  // Test insertion of entries
  check(n_ranks, true);
  check(2 * n_ranks, true);

  // Test empty sparsity pattern
  check(n_ranks, false);

  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    deallog << "OK" << std::endl;
}
