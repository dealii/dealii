// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2017 - 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

#include <deal.II/base/exception_macros.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/types.h>

#include <deal.II/lac/psblas_vector.h>

#include <cmath>
#include <vector>

#include "../../tests/tests.h"

using namespace dealii;

// The method reinit() rebuilds the PSBLAS descriptor on *all* processes when
// the partitioning changes on only some of them.
// The global size stays at 10 while the local sizes go from (2, 4, 4) to (2,
// 5, 3), so process 0 sees no change at all. Allocating a descriptor is a
// collective operation: if the processes decide on their own whether to do so,
// this test would deadlocks.

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  const MPI_Comm                   mpi_communicator = MPI_COMM_WORLD;
  AssertThrow(Utilities::MPI::n_mpi_processes(mpi_communicator) == 3,
              ExcMessage("This test needs to be run with 3 MPI processes."));

  MPILogInitAll log;

  const unsigned int rank = Utilities::MPI::this_mpi_process(mpi_communicator);

  const types::global_dof_index n_dofs = 10;

  const types::global_dof_index first_begin[]  = {0, 2, 6};
  const types::global_dof_index first_end[]    = {2, 6, 10};
  const types::global_dof_index second_begin[] = {0, 2, 7};
  const types::global_dof_index second_end[]   = {2, 7, 10};

  IndexSet first_partitioning(n_dofs);
  first_partitioning.add_range(first_begin[rank], first_end[rank]);
  IndexSet second_partitioning(n_dofs);
  second_partitioning.add_range(second_begin[rank], second_end[rank]);

  // Every process also reads the first entry owned by its right neighbor. For
  // process 0 that index stays the same, so its halo does not change either.
  IndexSet first_relevant(first_partitioning);
  IndexSet second_relevant(second_partitioning);
  if (rank + 1 < Utilities::MPI::n_mpi_processes(mpi_communicator))
    {
      first_relevant.add_index(first_end[rank]);
      second_relevant.add_index(second_end[rank]);
    }

  // Each entry is set to its global index, so the sum over all of them is
  // independent of the partitioning.
  const auto fill = [](PSCToolkitWrappers::Vector &v,
                       const IndexSet             &partitioning) {
    v = 0.;
    for (const types::global_dof_index i : partitioning)
      v(i) += static_cast<double>(i);
    v.compress(VectorOperation::add);
  };

  const auto check = [&](const PSCToolkitWrappers::Vector &owned,
                         const PSCToolkitWrappers::Vector &ghosted,
                         const types::global_dof_index     ghost_index) {
    AssertThrow(owned.size() == n_dofs,
                ExcMessage("The vector has the wrong global size."));
    AssertThrow(std::abs(owned.l1_norm() - 45.) < 1.e-12,
                ExcMessage("The vector does not contain the expected values."));

    if (ghost_index < n_dofs)
      {
        const std::vector<types::global_dof_index> indices = {ghost_index};
        std::vector<double>                        values(1);
        ghosted.extract_subvector_to(indices.begin(),
                                     indices.end(),
                                     values.begin());
        AssertThrow(std::abs(values[0] - static_cast<double>(ghost_index)) <
                      1.e-12,
                    ExcMessage("The ghost entry has the wrong value."));
      }

    deallog << "OK" << std::endl;
  };

  PSCToolkitWrappers::Vector owned(first_partitioning, mpi_communicator);
  PSCToolkitWrappers::Vector ghosted(first_partitioning,
                                     first_relevant,
                                     mpi_communicator);

  fill(owned, first_partitioning);
  ghosted = owned;
  check(owned, ghosted, first_end[rank]);

  // Repartition. Process 0 keeps both its owned and its ghost indices.
  owned.reinit(second_partitioning, mpi_communicator);
  ghosted.reinit(second_partitioning, second_relevant, mpi_communicator);

  fill(owned, second_partitioning);
  ghosted = owned;
  check(owned, ghosted, second_end[rank]);

  // ... and back to the original partitioning.
  owned.reinit(first_partitioning, mpi_communicator);
  ghosted.reinit(first_partitioning, first_relevant, mpi_communicator);

  fill(owned, first_partitioning);
  ghosted = owned;
  check(owned, ghosted, first_end[rank]);

  return 0;
}
