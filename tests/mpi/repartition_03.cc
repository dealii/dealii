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



// test that a second call to repartition() doesn't do anything

#include "../tests.h"
#include <deal.II/base/tensor.h>
#include <deal.II/grid/tria.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/base/utilities.h>



template <int dim>
void
print_cells(parallel::distributed::Triangulation<dim> &tr)
{
  for (typename Triangulation<dim>::active_cell_iterator
       cell = tr.begin_active();
       cell != tr.end(); ++cell)
    if (cell->is_locally_owned())
      deallog << cell->id() << std::endl;
}

template <int dim>
void
test()
{
  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);

  if (true)
    {
      parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD,
                                                   dealii::Triangulation<dim,dim>::none,
                                                   parallel::distributed::Triangulation<dim>::no_automatic_repartitioning);

      GridGenerator::hyper_cube(tr);
      tr.refine_global(2);

      deallog << "*** 1. everything on one core:" << std::endl;

      deallog << "locally owned cells: " << tr.n_locally_owned_active_cells()
              << " / "
              << tr.n_global_active_cells()
              << std::endl;
      print_cells(tr);

      deallog << "*** 2. repartition:" << std::endl;

      tr.repartition();

      deallog << "locally owned cells: " << tr.n_locally_owned_active_cells()
              << " / "
              << tr.n_global_active_cells()
              << std::endl;

      print_cells(tr);

      deallog << "*** 3. repartition again (noop):" << std::endl;
      tr.repartition();

      deallog << "locally owned cells: " << tr.n_locally_owned_active_cells()
              << " / "
              << tr.n_global_active_cells()
              << std::endl;

      print_cells(tr);

      const unsigned int checksum = tr.get_checksum ();
      if (myid == 0)
        deallog << "Checksum: "
                << checksum
                << std::endl;
    }

  if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
    deallog << "OK" << std::endl;
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv, 1);
  MPILogInitAll log;
  test<2>();
}
