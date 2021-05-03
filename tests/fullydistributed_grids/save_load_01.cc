// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2021 by the deal.II authors
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



// Test fullydistributed::Triangulation::load()/save().
// Create a triangulation, save it and load it.
// The initial and the loaded triangulations must be the same.

#include <deal.II/distributed/fully_distributed_tria.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_description.h>

#include <deal.II/lac/la_parallel_vector.h>

#include "./tests.h"

using namespace dealii;

template <int dim, typename TriangulationType>
void
test(TriangulationType &triangulation)
{
  const std::string filename = "save_load_" + std::to_string(dim) + "d_out";

  triangulation.save(filename);
  deallog.push("save");
  print_statistics(triangulation, false);
  deallog.pop();

  triangulation.clear();

  triangulation.load(filename);
  deallog.push("load");
  print_statistics(triangulation, false);
  deallog.pop();
}


int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  MPILogInitAll all;

  deallog.push("2d");
  {
    constexpr int dim = 2;

    parallel::distributed::Triangulation<dim> triangulation(MPI_COMM_WORLD);
    GridGenerator::hyper_cube(triangulation);
    triangulation.refine_global(3);

    const auto description = TriangulationDescription::Utilities::
      create_description_from_triangulation(triangulation, MPI_COMM_WORLD);

    parallel::fullydistributed::Triangulation<dim> triangulation_pft(
      MPI_COMM_WORLD);
    triangulation_pft.create_triangulation(description);

    test<dim>(triangulation_pft);
  }
  deallog.pop();

  deallog.push("3d");
  {
    constexpr int dim = 3;

    parallel::distributed::Triangulation<dim> triangulation(MPI_COMM_WORLD);
    GridGenerator::hyper_cube(triangulation);
    triangulation.refine_global(3);

    const auto description = TriangulationDescription::Utilities::
      create_description_from_triangulation(triangulation, MPI_COMM_WORLD);

    parallel::fullydistributed::Triangulation<dim> triangulation_pft(
      MPI_COMM_WORLD);
    triangulation_pft.create_triangulation(description);

    test<dim>(triangulation_pft);
  }
  deallog.pop();
}
