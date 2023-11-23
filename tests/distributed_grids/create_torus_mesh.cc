// ---------------------------------------------------------------------
//
// Copyright (C) 2023 - 2023 by the deal.II authors
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



// Test GridGenerator::torus, which used to run into an assertion with p4est
// due to bad orientations before 2023, see
// https://github.com/dealii/dealii/issues/16272

#include <deal.II/distributed/tria.h>

#include <deal.II/grid/grid_generator.h>

#include "../tests.h"


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  initlog();

  // create mesh
  const int dim = 3;

  parallel::distributed::Triangulation<dim> triangulation(
    MPI_COMM_WORLD,
    typename Triangulation<dim>::MeshSmoothing(
      Triangulation<dim>::smoothing_on_refinement |
      Triangulation<dim>::smoothing_on_coarsening));

  const double       R                = 2.;
  const double       r                = 1.;
  const unsigned int n_cells_toroidal = 1;
  const double       phi              = 0.5 * dealii::numbers::PI;
  dealii::GridGenerator::torus<dim, dim>(
    triangulation, R, r, n_cells_toroidal, phi);
  triangulation.reset_all_manifolds();
  triangulation.refine_global(1);

  deallog << "Tria info: " << triangulation.n_vertices() << " "
          << triangulation.n_cells() << " " << triangulation.n_faces() << " "
          << triangulation.n_raw_lines() << std::endl;

  return 0;
}
