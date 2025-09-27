// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// DoFTools::locally_owned_dofs_per_subdomain wrongly sized the array
// it returns when we have some processors that don't own any cells

#include <deal.II/distributed/shared_tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>

#include "../tests.h"


template <int dim>
void
test()
{
  parallel::shared::Triangulation<dim> triangulation(
    MPI_COMM_WORLD,
    ::Triangulation<dim>::none,
    false,
    parallel::shared::Triangulation<dim>::partition_zorder);
  GridGenerator::hyper_cube(triangulation, -1.0, 1.0);

  FE_Q<dim>       fe(1);
  DoFHandler<dim> dof_handler(triangulation);
  dof_handler.distribute_dofs(fe);

  // this used to crash here:
  const IndexSet locally_owned_dofs = dof_handler.locally_owned_dofs();

  deallog << "dim=" << dim << std::endl
          << "rank=" << Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)
          << std::endl
          << "  n_dofs=" << dof_handler.n_dofs() << std::endl
          << "  n_locally_owned_dofs=" << locally_owned_dofs.n_elements()
          << std::endl;
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());

  MPILogInitAll all;

  test<2>();
  test<3>();
}
