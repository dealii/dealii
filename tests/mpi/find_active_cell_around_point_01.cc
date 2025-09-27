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



// make sure only one processor finds a locally-owned cell around a point

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
      if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
        deallog << "hyper_cube" << std::endl;

      parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);

      GridGenerator::hyper_cube(tr);
      tr.refine_global(2);

      // choose a point that is guaranteed to lie in the domain but not
      // at the interface between cells
      Point<dim> p;
      for (unsigned int d = 0; d < dim; ++d)
        p[d] = 1. / 3;

      typename parallel::distributed::Triangulation<dim>::active_cell_iterator
        cell = GridTools::find_active_cell_around_point(tr, p);

      const unsigned int n_locally_owned =
        Utilities::MPI::sum(cell->is_locally_owned() ? 1 : 0, MPI_COMM_WORLD);

      const unsigned int n_locally_owned_or_ghost =
        Utilities::MPI::sum(!cell->is_artificial() ? 1 : 0, MPI_COMM_WORLD);

      if (myid == 0)
        deallog << "Locally owned: " << n_locally_owned << std::endl
                << "Locally owned or ghost: " << n_locally_owned_or_ghost
                << std::endl;
    }
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);


  deallog.push(Utilities::int_to_string(myid));

  if (myid == 0)
    {
      initlog();

      deallog.push("2d");
      test<2>();
      deallog.pop();
      deallog.push("3d");
      test<3>();
      deallog.pop();
    }
  else
    {
      test<2>();
      test<3>();
    }
}
