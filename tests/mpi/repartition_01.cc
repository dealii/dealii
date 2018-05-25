// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2017 by the deal.II authors
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



// test manual repartitioning

#include <deal.II/base/tensor.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>

#include "../tests.h"



template <int dim>
void
test()
{
  unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  if (true)
    {
      parallel::distributed::Triangulation<dim> tr(
        MPI_COMM_WORLD,
        dealii::Triangulation<dim, dim>::none,
        parallel::distributed::Triangulation<dim>::no_automatic_repartitioning);

      GridGenerator::hyper_cube(tr);
      tr.refine_global(2);

      deallog << "locally owned cells: " << tr.n_locally_owned_active_cells()
              << " / " << tr.n_global_active_cells() << std::endl;

      // tr.write_mesh_vtk("a");

      tr.repartition();

      // tr.write_mesh_vtk("b");

      deallog << "locally owned cells: " << tr.n_locally_owned_active_cells()
              << " / " << tr.n_global_active_cells() << std::endl;

      const unsigned int checksum = tr.get_checksum();
      if (myid == 0)
        deallog << "Checksum: " << checksum << std::endl;
    }

  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    deallog << "OK" << std::endl;
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;
  test<2>();
}
