// ---------------------------------------------------------------------
//
// Copyright (C) 2021 by the deal.II authors
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


// Create a serial triangulation with empty cell_info fields and copy it.

#include <deal.II/base/mpi.h>

#include <deal.II/distributed/fully_distributed_tria.h>
#include <deal.II/distributed/shared_tria.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_description.h>

#include "./tests.h"

using namespace dealii;

template <int dim>
void
test(MPI_Comm comm)
{
  // create serial triangulation
  Triangulation<dim> basetria;
  GridGenerator::hyper_cube(basetria);

  // create instance of pft
  parallel::fullydistributed::Triangulation<dim> tria_pft(comm);

  // extract relevant information form serial triangulation
  auto construction_data =
    TriangulationDescription::Utilities::create_description_from_triangulation(
      basetria, comm);

  // create some empty levels to imitate empty levels on some MPI processes
  construction_data.cell_infos.resize(basetria.n_levels() + 2);

  // actually create triangulation
  tria_pft.create_triangulation(construction_data);

  // if everything went fine, we do not have bugs inside the library
  deallog << "OK" << std::endl;
}

int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  mpi_initlog();

  const MPI_Comm comm = MPI_COMM_WORLD;

  test<2>(comm);
}
