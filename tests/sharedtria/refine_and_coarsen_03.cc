// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------



// Check that calling prepare_coarsening_and_refinement() directly on a
// parallel::shared::Triangulation, as is customary before attaching data
// for a solution transfer, synchronizes refinement and coarsening flags
// over all processes *before* flag smoothing runs.
//
// Refinement and coarsening flags only need to be set on locally owned
// cells. If prepare_coarsening_and_refinement() does not synchronize the
// flags first, smoothing (in particular the rule that a family of sibling
// cells is only coarsened if all children are flagged) only sees the
// flags this process has set itself. Flags on cells owned by other
// processes, which are artificial here, default to false. For a family
// that straddles a partition boundary this drops the coarsen flag on
// every process, and the resulting mesh ends up depending on the
// partitioning.
//
// We check this by counting, immediately after
// prepare_coarsening_and_refinement(), how many active cells (including
// artificial ones) have a refine/coarsen flag set: with synchronized
// flags this count has to be identical on every process. The cell counts
// logged by this test therefore have to agree between the mpirun=1 and
// mpirun=3 output files.

#include <deal.II/distributed/shared_tria.h>

#include <deal.II/grid/grid_generator.h>

#include "../tests.h"



template <int dim>
void
test()
{
  parallel::shared::Triangulation<dim> tria(
    MPI_COMM_WORLD,
    Triangulation<dim>::limit_level_difference_at_vertices,
    /*allow_artificial_cells*/ true,
    parallel::shared::Triangulation<dim>::partition_zorder);

  GridGenerator::subdivided_hyper_cube(tria, 2);
  tria.refine_global(2);
  deallog << "cells before: " << tria.n_global_active_cells() << std::endl;

  // Set refinement and coarsening flags on locally owned cells only:

  for (const auto &cell : tria.active_cell_iterators())
    if (cell->is_locally_owned())
      {
        const auto center = cell->center();

        bool refine = true;
        for (unsigned int d = 0; d < dim; ++d)
          refine = refine && (center[d] < 0.25);

        if (refine)
          cell->set_refine_flag();
        else if (center[0] > 0.5)
          cell->set_coarsen_flag();
      }

  tria.prepare_coarsening_and_refinement();

  // Immediately after prepare_coarsening_and_refinement() every process
  // must see the complete, partition-independent set of flags. We count
  // over *all* active cells (including artificial ones), since a flag on
  // an artificial cell can only ever become visible through
  // synchronization:

  unsigned int n_refine_flags  = 0;
  unsigned int n_coarsen_flags = 0;
  for (const auto &cell : tria.active_cell_iterators())
    {
      if (cell->refine_flag_set())
        ++n_refine_flags;
      if (cell->coarsen_flag_set())
        ++n_coarsen_flags;
    }
  deallog << "flags after prepare: refine=" << n_refine_flags
          << " coarsen=" << n_coarsen_flags << std::endl;

  AssertThrow(Utilities::MPI::max(n_refine_flags, MPI_COMM_WORLD) ==
                Utilities::MPI::min(n_refine_flags, MPI_COMM_WORLD),
              ExcInternalError());
  AssertThrow(Utilities::MPI::max(n_coarsen_flags, MPI_COMM_WORLD) ==
                Utilities::MPI::min(n_coarsen_flags, MPI_COMM_WORLD),
              ExcInternalError());

  tria.execute_coarsening_and_refinement();
  deallog << "cells after: " << tria.n_global_active_cells() << std::endl;

  deallog << "OK" << std::endl;
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  mpi_initlog();

  deallog.push("1d");
  test<1>();
  deallog.pop();

  deallog.push("2d");
  test<2>();
  deallog.pop();

  deallog.push("3d");
  test<3>();
  deallog.pop();
}
