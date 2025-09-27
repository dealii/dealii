// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


/**
 * Test fix for bug (no ghost value update in the case of DG)
 * exposed in https://github.com/dealii/dealii/pull/16049.
 */

#include <deal.II/base/partitioner.h>

#include <deal.II/distributed/shared_tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_dgq.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/multigrid/mg_transfer_global_coarsening.h>

#include "../tests.h"

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    all;

  const int dim = 2;
  using Number  = double;

  // fine mesh
  parallel::shared::Triangulation<dim> tria_fine(
    MPI_COMM_WORLD,
    ::Triangulation<dim>::none,
    true,
    parallel::shared::Triangulation<dim>::partition_custom_signal);
  tria_fine.signals.create.connect([&]() {
    for (const auto &cell : tria_fine.active_cell_iterators())
      cell->set_subdomain_id(cell->active_cell_index());
  });
  GridGenerator::subdivided_hyper_rectangle(tria_fine,
                                            {2, 1},
                                            Point<dim>(0, 0),
                                            Point<dim>(2, 1));

  DoFHandler<dim> dof_handler_fine(tria_fine);
  dof_handler_fine.distribute_dofs(FE_DGQ<dim>(1));

  // coarse mesh
  parallel::shared::Triangulation<dim> tria_coarse(
    MPI_COMM_WORLD,
    ::Triangulation<dim>::none,
    true,
    parallel::shared::Triangulation<dim>::partition_custom_signal);
  tria_coarse.signals.create.connect([&]() {
    for (const auto &cell : tria_coarse.active_cell_iterators())
      cell->set_subdomain_id(1 - cell->active_cell_index());
  });
  GridGenerator::subdivided_hyper_rectangle(tria_coarse,
                                            {2, 1},
                                            Point<dim>(0, 0),
                                            Point<dim>(2, 1));

  DoFHandler<dim> dof_handler_coarse(tria_coarse);
  dof_handler_coarse.distribute_dofs(FE_DGQ<dim>(1));

  // set up two-level transfer
  MGTwoLevelTransfer<dim, LinearAlgebra::distributed::Vector<Number>> transfer;
  transfer.reinit(dof_handler_fine, dof_handler_coarse);

  IndexSet is_all(dof_handler_coarse.n_dofs());
  is_all.add_range(0, dof_handler_coarse.n_dofs());

  const auto partitioner_coarse = std::make_shared<Utilities::MPI::Partitioner>(
    dof_handler_coarse.locally_owned_dofs(), is_all, MPI_COMM_WORLD);

  const auto partitioner_fine = std::make_shared<Utilities::MPI::Partitioner>(
    dof_handler_fine.locally_owned_dofs(), is_all, MPI_COMM_WORLD);

  transfer.enable_inplace_operations_if_possible(partitioner_coarse,
                                                 partitioner_fine);

  LinearAlgebra::distributed::Vector<Number> vec_fine(partitioner_fine);
  LinearAlgebra::distributed::Vector<Number> vec_coarse(partitioner_coarse);

  transfer.restrict_and_add(vec_coarse, vec_fine);

  deallog << "OK!" << std::endl;
}
