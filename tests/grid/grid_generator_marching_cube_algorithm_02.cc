// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test GridTools::MarchingCubeAlgorithm<dim>::process()

#include <deal.II/base/function_signed_distance.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

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


using VectorType = LinearAlgebra::distributed::Vector<double>;

template <int dim>
class HeavisideFunction : public Function<dim>
{
public:
  HeavisideFunction(const Point<dim> &center, const double radius)
    : Function<dim>(1)
    , distance_sphere(center, radius)
  {}

  double
  value(const Point<dim> &p, const unsigned int component) const
  {
    const double distance = distance_sphere.value(p);

    // compute sign
    return boost::math::sign(distance);
  }

  Functions::SignedDistance::Sphere<dim> distance_sphere;
};

template <int dim>
void
create_mca_tria(const unsigned int   n_subdivisions,
                const double         iso_level,
                const Function<dim> &my_function)
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

  VectorType     ls_vector;
  const IndexSet locally_relevant_dofs =
    DoFTools::extract_locally_relevant_dofs(background_dof_handler);


  ls_vector.reinit(background_dof_handler.locally_owned_dofs(),
                   locally_relevant_dofs,
                   MPI_COMM_WORLD);

  dealii::VectorTools::interpolate(mapping,
                                   background_dof_handler,
                                   my_function,
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
                                          triangulation.get_mpi_communicator());
    }
}

template <int dim>
void
test()
{
  for (unsigned int i = 1; i <= 3; ++i)
    {
      {
        deallog << "signed distance function" << std::endl;
        const auto my_function =
          Functions::SignedDistance::Sphere<dim>(Point<dim>(), 0.75);
        create_mca_tria<dim>(i, -0.1 + i * 0.05, my_function);
      }
      {
        // test function with constant regions
        deallog << "heaviside function" << std::endl;
        const auto my_function = HeavisideFunction<dim>(Point<dim>(), 0.75);
        create_mca_tria<dim>(i, 0, my_function);
      }
    }
}

int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    all;

  deallog << std::scientific << std::setprecision(6);

  test<1>();
  test<2>();

  return 0;
}
