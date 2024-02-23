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



// Test interaction with p4est with a complicated 3d grid read from file.
//
// the files we read here are the ones already used in deal.II/grid_in_3d

#include <deal.II/base/tensor.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"

#include "coarse_grid_common.h"



template <int dim>
void
test(const char *filename)
{
  const char *p = strrchr(filename, '/');
  deallog.push(p);

  parallel::distributed::Triangulation<dim> tr(
    MPI_COMM_WORLD,
    Triangulation<dim>::none,
    parallel::distributed::Triangulation<dim>::communicate_vertices_to_p4est);

  GridIn<dim> gi;
  gi.attach_triangulation(tr);
  std::ifstream in(filename);
  gi.read_xda(in);

  write_vtk(tr, "1");

  deallog.pop();
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  initlog();

  deallog.push("3d");

  test<3>(SOURCE_DIR "/../grid/grid_in_3d/1.in");
  test<3>(SOURCE_DIR "/../grid/grid_in_3d/2.in");
  test<3>(SOURCE_DIR "/../grid/grid_in_3d/3.in");
  test<3>(SOURCE_DIR "/../grid/grid_in_3d/4.in");

  test<3>(SOURCE_DIR "/../grid/grid_in_3d/evil_0.in");
  test<3>(SOURCE_DIR "/../grid/grid_in_3d/evil_1.in");
  test<3>(SOURCE_DIR "/../grid/grid_in_3d/evil_2.in");
  test<3>(SOURCE_DIR "/../grid/grid_in_3d/evil_3.in");
  test<3>(SOURCE_DIR "/../grid/grid_in_3d/evil_4.in");

  deallog.pop();
}
