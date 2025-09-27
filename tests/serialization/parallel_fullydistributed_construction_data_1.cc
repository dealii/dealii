// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// check serialization for TriangulationDescription::Description<dim, spacedim>

#include <deal.II/distributed/fully_distributed_tria.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_description.h>

#include <boost/serialization/vector.hpp>

#include "../grid/tests.h"
#include "serialization.h"


template <int dim>
void
test(MPI_Comm comm)
{
  // create serial triangulation
  Triangulation<dim> basetria;
  GridGenerator::hyper_cube(basetria);
  basetria.refine_global(1);

  GridTools::partition_triangulation_zorder(
    Utilities::MPI::n_mpi_processes(comm), basetria);

  auto t1 =
    TriangulationDescription::Utilities::create_description_from_triangulation(
      basetria, comm);

  // compare equal TriangulationDescription::Descriptions
  auto t2 = t1;
  verify(t1, t2);

  basetria.refine_global(1);

  GridTools::partition_triangulation_zorder(
    Utilities::MPI::n_mpi_processes(comm), basetria);

  auto t3 =
    TriangulationDescription::Utilities::create_description_from_triangulation(
      basetria, comm);

  // compare different TriangulationDescription::Descriptions
  verify(t1, t3);
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    all;
  deallog << std::setprecision(3);
  const MPI_Comm comm = MPI_COMM_WORLD;

  test<1>(comm);
  test<2>(comm);
  test<3>(comm);

  deallog << "OK" << std::endl;
}
