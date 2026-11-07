// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2023 - 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


/**
 * Test the class Functions::SignedDistance::ArbitraryLevelSet by creating an
 * object and evaluating the level set at points where the expected values
 * are known.
 */

#include <deal.II/base/function_signed_distance.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"

// set this flag to true to enable vtu output for visualizing the level set.
constexpr bool enable_vtu = false;
namespace
{
  template <int dim>
  void
  print_values_at_point(const dealii::Function<dim> &function,
                        const dealii::Point<dim>    &point)
  {
    deallog << "point = " << point << std::endl;
    deallog << "value " << function.value(point) << std::endl;
  }

  template <int dim>
  void
  output_signed_distance_vtu(const dealii::Function<dim> &signed_distance,
                             const std::string           &filename)
  {
    if (not enable_vtu)
      return;

    dealii::Triangulation<dim> triangulation;
    dealii::GridGenerator::hyper_cube(triangulation, -3, 3);
    triangulation.refine_global(5);

    dealii::FE_Q<dim>       fe(1);
    dealii::DoFHandler<dim> dof_handler(triangulation);
    dof_handler.distribute_dofs(fe);

    dealii::Vector<double> values(dof_handler.n_dofs());

    dealii::VectorTools::interpolate(dof_handler, signed_distance, values);

    dealii::DataOut<dim> data_out;
    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(values, "signed_distance");

    data_out.build_patches();

    std::ofstream output(filename);
    data_out.write_vtu(output);
  }

  template <int dim>
  void
  test_arbitrary_level_set_signed_distance_sphere()
  {
    deallog << "test_arbitrary_level_set_sphere" << std::endl;

    // the images are binary (1 if inside, 0 if outside)
    // choose value larger than 0 to ignore the void
    const float min_density = 0.1;
    // choose value larger than 1 to capture the sphere
    const float max_density = 1.1;

    // create the level set.
    // It's a sphere in [0,100]^3 with a center at (50, 50, 50) and radius
    // of 25.
    const Functions::SignedDistance::ArbitraryLevelSet<dim>
      signed_distance_sphere(std::string(SOURCE_DIR
                                         "/arbitrary_level_set/sphere_3d.nii"),
                             min_density,
                             max_density);

    dealii::Point<dim> p(50.0, 50.0, 50.0);


    deallog << "center" << std::endl;
    print_values_at_point(signed_distance_sphere, p);
    p[0] = 0;
    print_values_at_point(signed_distance_sphere, p);
    p[0] = 25;
    print_values_at_point(signed_distance_sphere, p);
    p[0] = 35;
    p[1] = 35;
    print_values_at_point(signed_distance_sphere, p);
  }

  template <int dim>
  void
  test_arbitrary_level_set_signed_distance_cube()
  {
    deallog << "test_arbitrary_level_set_cube" << std::endl;

    // the images are binary (1 if inside, 0 if outside)
    // choose value larger than 0 to ignore the void
    const float min_density = 0.1;
    // choose value larger than 1 to capture the cube
    const float max_density = 1.1;

    // create the level set.
    // It's a cube in [0,100]^3 with a center at (50, 50, 50) and a side length
    // of 50.
    const dealii::Functions::SignedDistance::ArbitraryLevelSet<dim>
      signed_distance_cube(std::string(SOURCE_DIR
                                       "/arbitrary_level_set/cube.nii"),
                           min_density,
                           max_density);

    dealii::Point<dim> p(50.0, 50.0, 50.0);


    deallog << "center" << std::endl;
    print_values_at_point(signed_distance_cube, p);
    p[0] = 0;
    print_values_at_point(signed_distance_cube, p);
    p[0] = 25;
    print_values_at_point(signed_distance_cube, p);
    p[0] = 35;
    p[1] = 35;
    print_values_at_point(signed_distance_cube, p);
  }


  template <int dim>
  void
  run_test()
  {
    deallog << "dim = " << dim << std::endl;
    deallog << std::endl;

    test_arbitrary_level_set_signed_distance_sphere<dim>();
    test_arbitrary_level_set_signed_distance_cube<dim>();

    deallog << std::endl;
  }
} // namespace

int
main()
{
  run_test<3>();
}
