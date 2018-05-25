// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2017 by the deal.II Authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

// test DoFRenumbering::block_wise with parallel::shared::Triangulation

#include "../tests.h"

// all include files you need here
#include <deal.II/distributed/shared_tria.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>

template <int dim>
void
test()
{
  const unsigned int this_mpi_process =
    Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  FESystem<dim> fe(FE_Q<dim>(1), 2);

  {
    parallel::shared::Triangulation<dim> tria(
      MPI_COMM_WORLD,
      Triangulation<dim>::none,
      false,
      parallel::shared::Triangulation<dim>::Settings::partition_zorder);
    GridGenerator::hyper_cube(tria, -1, 1);
    tria.refine_global(2);

    DoFHandler<dim> dh(tria);
    dh.distribute_dofs(fe);

    DoFRenumbering::block_wise(dh);

    const std::vector<IndexSet> locally_owned_dofs_per_subdomain =
      DoFTools::locally_owned_dofs_per_subdomain(dh);

    const types::global_dof_index dofs_per_block = dh.n_dofs() / 2;
    std::vector<IndexSet>         locally_owned_dofs_per_block(2);
    locally_owned_dofs_per_block[0] =
      locally_owned_dofs_per_subdomain[this_mpi_process].get_view(
        0, dofs_per_block);
    locally_owned_dofs_per_block[1] =
      locally_owned_dofs_per_subdomain[this_mpi_process].get_view(
        dofs_per_block, dh.n_dofs());

    if (locally_owned_dofs_per_block[0] != locally_owned_dofs_per_block[1])
      AssertThrow(false,
                  ExcMessage("Locally owned dofs differ across blocks."));
  }

  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    deallog << "OK for " << dim << "d" << std::endl;
}

int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;

  test<2>();
  test<3>();

  return 0;
}
