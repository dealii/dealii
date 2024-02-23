// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// load a triangulation from an existing forest


#include <deal.II/distributed/tria.h>

#include <deal.II/grid/grid_generator.h>

#include "../tests.h"


template <int dim>
void
test()
{
  // prepare two independent triangulations with the same coarse mesh
  parallel::distributed::Triangulation<dim> tria1(MPI_COMM_WORLD),
    tria2(MPI_COMM_WORLD);
  GridGenerator::subdivided_hyper_cube(tria1, 2);
  GridGenerator::subdivided_hyper_cube(tria2, 2);

  // refine the first, and transfer the forest to the other
  tria1.refine_global(2);
  tria2.load(tria1.get_p4est());

  // verify that global information matches
  AssertDimension(tria1.n_global_levels(), tria2.n_global_levels());
  AssertDimension(tria1.n_global_active_cells(), tria2.n_global_active_cells());

  deallog << "OK" << std::endl;
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;

  deallog.push("2d");
  test<2>();
  deallog.pop();
  deallog.push("3d");
  test<3>();
  deallog.pop();
}
