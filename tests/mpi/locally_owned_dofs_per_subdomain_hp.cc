// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------


// DoFTools::locally_owned_dofs_per_subdomain wrongly sized the array
// it returns when we have some processors that don't own any cells
//
// like the test without _hp, but for a hp::DoFHandler

#include <deal.II/distributed/shared_tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/hp/fe_collection.h>

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

  hp::FECollection<dim> fe;
  fe.push_back(FE_Q<dim>(1));
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
