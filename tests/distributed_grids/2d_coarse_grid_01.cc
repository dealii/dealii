// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Test interaction with p4est with a few simple coarse grids in 2d

#include <deal.II/base/tensor.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"

#include "coarse_grid_common.h"



template <int dim>
void
test(std::ostream & /*out*/)
{
  if (true)
    {
      deallog << "hyper_cube" << std::endl;

      parallel::distributed::Triangulation<dim> tr(
        MPI_COMM_WORLD,
        Triangulation<dim>::none,
        parallel::distributed::Triangulation<
          dim>::communicate_vertices_to_p4est);

      GridGenerator::hyper_cube(tr);
      write_vtk(tr, "1");
    }


  if (true)
    {
      deallog << "hyper_ball" << std::endl;

      parallel::distributed::Triangulation<dim> tr(
        MPI_COMM_WORLD,
        Triangulation<dim>::none,
        parallel::distributed::Triangulation<
          dim>::communicate_vertices_to_p4est);

      GridGenerator::hyper_ball(tr, Point<dim>(), 3.);
      write_vtk(tr, "2");
    }

  if (true)
    {
      deallog << "half_hyper_ball" << std::endl;

      parallel::distributed::Triangulation<dim> tr(
        MPI_COMM_WORLD,
        Triangulation<dim>::none,
        parallel::distributed::Triangulation<
          dim>::communicate_vertices_to_p4est);

      GridGenerator::half_hyper_ball(tr, Point<dim>(), 3.);
      write_vtk(tr, "3");
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
