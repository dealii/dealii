// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2020 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


// Test that
// MGTransferGlobalCoarseningTools::create_geometric_coarsening_sequence works
// in case the initial parallel distributed triangulation disables automatic
// repartitioning.

#include <deal.II/base/mpi.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/multigrid/mg_transfer_global_coarsening.h>

#include "../tests.h"

template <int dim>
void
test()
{
  parallel::distributed::Triangulation<dim> tria(
    MPI_COMM_WORLD,
    Triangulation<dim>::limit_level_difference_at_vertices,
    parallel::distributed::Triangulation<dim>::no_automatic_repartitioning);

  GridGenerator::hyper_cube(tria);
  tria.refine_global(3);
  tria.repartition();

  deallog << "Test in " << dim << "d" << std::endl;
  auto triangulations =
    MGTransferGlobalCoarseningTools::create_geometric_coarsening_sequence(tria);
  for (const auto tr : triangulations)
    deallog << "Active cells in level triangulation: "
            << tr->n_global_active_cells() << std::endl;
  deallog << std::endl;
}

int
main(int argc, char **argv)
{
  mpi_initlog();
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

  test<2>();
  test<3>();
}
