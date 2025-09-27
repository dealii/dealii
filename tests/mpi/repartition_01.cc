// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



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
