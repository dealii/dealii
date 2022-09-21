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

// Test GridTools::MarchingCubeAlgorithm<dim>::process()

#include <deal.II/base/function_signed_distance.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/lac/la_parallel_vector.h>

#include <deal.II/matrix_free/matrix_free.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <iostream>

#include "../tests.h"

using namespace dealii;

using VectorType = LinearAlgebra::distributed::Vector<double>;

template <int dim>
void
test(const unsigned int n_subdivisions, const double iso_level)
{
  deallog << "dim=" << dim << " iso level: " << iso_level << std::endl;
  const int degree        = 3;
  const int n_refinements = 5;

  parallel::shared::Triangulation<dim> triangulation(MPI_COMM_WORLD);

  GridGenerator::hyper_cube(triangulation, -1.5, 1.5);

  triangulation.refine_global(n_refinements);

  FE_Q<dim>      fe(degree);
  MappingQ1<dim> mapping;

  DoFHandler<dim> background_dof_handler;
  background_dof_handler.reinit(triangulation);
  background_dof_handler.distribute_dofs(fe);

  VectorType ls_vector;
  IndexSet   locally_relevant_dofs;
  DoFTools::extract_locally_relevant_dofs(background_dof_handler,
                                          locally_relevant_dofs);


  ls_vector.reinit(background_dof_handler.locally_owned_dofs(),
                   locally_relevant_dofs,
                   MPI_COMM_WORLD);

  dealii::VectorTools::interpolate(
    mapping,
    background_dof_handler,
    Functions::SignedDistance::Sphere<dim>(Point<dim>(), 0.75),
    ls_vector);

  std::vector<Point<dim>> vertices;

  const GridTools::MarchingCubeAlgorithm<dim, VectorType> mc(
    mapping, background_dof_handler.get_fe(), n_subdivisions);

  ls_vector.update_ghost_values();
  mc.process(background_dof_handler, ls_vector, iso_level, vertices);

  for (const auto &v : vertices)
    deallog << "point found: " << v << std::endl;

  if (false /*write background mesh*/)
    {
      DataOut<dim> data_out;
      data_out.attach_dof_handler(background_dof_handler);
      data_out.add_data_vector(ls_vector, "level_set");
      data_out.build_patches(2);
      data_out.write_vtu_with_pvtu_record("./",
                                          "data_background_" +
                                            std::to_string(n_subdivisions),
                                          0,
                                          triangulation.get_communicator());
    }
}

int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    all;

  deallog << std::scientific << std::setprecision(6);

  // dim==1
  for (unsigned int i = 1; i <= 3; ++i)
    test<1>(i, -0.1 + i * 0.05);

  // dim==2
  for (unsigned int i = 1; i <= 3; ++i)
    test<2>(i, -0.1 + i * 0.05);

  return 0;
}
