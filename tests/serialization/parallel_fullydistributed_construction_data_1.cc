// ---------------------------------------------------------------------
//
// Copyright (C) 2019 by the deal.II authors
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


// check serialization for TriangulationDescription::Description<dim, spacedim>

#include <deal.II/distributed/fully_distributed_tria.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_description.h>

#include <boost/serialization/vector.hpp>

#include "../fullydistributed_grids/tests.h"
#include "serialization.h"

using namespace dealii;

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
