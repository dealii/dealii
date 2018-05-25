// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2017 by the deal.II authors
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


// create a shared tria mesh and copy it

#include <deal.II/distributed/shared_tria.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"

template <int dim>
void
test()
{
  Triangulation<dim> tr1;

  GridGenerator::hyper_cube(tr1);
  tr1.refine_global(2);

  deallog << " n_active_cells: " << tr1.n_active_cells() << "\n" << std::endl;

  parallel::shared::Triangulation<dim> tr2(
    MPI_COMM_WORLD,
    ::Triangulation<dim>::none,
    false,
    parallel::shared::Triangulation<dim>::partition_metis);
  tr2.copy_triangulation(tr1);

  deallog << " n_active_cells: " << tr2.n_active_cells() << "\n"
          << " locally_owned_subdomain(): " << tr2.locally_owned_subdomain()
          << "\n"
          << " n_locally_owned_active_cells: "
          << tr2.n_locally_owned_active_cells() << "\n"
          << " n_global_active_cells: " << tr2.n_global_active_cells() << "\n"
          << std::endl;

  parallel::shared::Triangulation<dim> tr3(
    MPI_COMM_WORLD,
    ::Triangulation<dim>::none,
    false,
    parallel::shared::Triangulation<dim>::partition_metis);
  tr3.copy_triangulation(tr2);

  deallog << " n_active_cells: " << tr3.n_active_cells() << "\n"
          << " locally_owned_subdomain(): " << tr3.locally_owned_subdomain()
          << "\n"
          << " n_locally_owned_active_cells: "
          << tr3.n_locally_owned_active_cells() << "\n"
          << " n_global_active_cells: " << tr3.n_global_active_cells() << "\n"
          << std::endl;
}

int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    all;

  deallog.push("2d");
  test<2>();
  deallog.pop();
  deallog.push("3d");
  test<3>();
  deallog.pop();
}
