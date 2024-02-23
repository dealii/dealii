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



// Like coarse_grid_x_01, but instead of checking that what we copy into
// the p4est data structure is correct, check that what we get back is
// correct

#include <deal.II/base/tensor.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"



template <int dim>
void
test(std::ostream & /*out*/)
{
  if (true)
    {
      deallog << "hyper_cube" << std::endl;

      parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);

      GridGenerator::hyper_cube(tr);

      deallog << "Triangulation copied back from p4est has "
              << tr.n_active_cells() << " active cells" << std::endl;
      GridOut().write_gnuplot(tr, deallog.get_file_stream());
    }


  if (true)
    {
      deallog << "hyper_ball" << std::endl;

      parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);

      GridGenerator::hyper_ball(tr, Point<dim>(), 3.);

      deallog << "Triangulation copied back from p4est has "
              << tr.n_active_cells() << " active cells" << std::endl;
      GridOut().write_gnuplot(tr, deallog.get_file_stream());
    }

  if (true)
    {
      deallog << "half_hyper_ball" << std::endl;

      parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);

      GridGenerator::half_hyper_ball(tr, Point<dim>(), 3.);

      deallog << "Triangulation copied back from p4est has "
              << tr.n_active_cells() << " active cells" << std::endl;
      GridOut().write_gnuplot(tr, deallog.get_file_stream());
    }
}


int
main(int argc, char *argv[])
{
  initlog();

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  deallog.push("2d");
  test<2>(deallog.get_file_stream());
  deallog.pop();
}
