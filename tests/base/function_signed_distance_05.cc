// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


/**
 * Test the class Functions::SignedDistance::InfiniteCylinder by creating an
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
  print_values_at_point(const Function<dim> &function, const Point<dim> &point)
  {
    deallog << "point = " << point << std::endl;
    deallog << "value = " << function.value(point) << std::endl;
  }



  template <int dim>
  void
  print_gradient_at_point(const Function<dim> &function,
                          const Point<dim>    &point)
  {
    deallog << "point = " << point << std::endl;
    deallog << "gradient = " << function.gradient(point) << std::endl;
  }



  template <int dim>
  void
  output_signed_distance_vtu(const Function<dim> &signed_distance,
                             const std::string   &filename)
  {
    if (not enable_vtu)
      return;

    Triangulation<dim> triangulation;
    GridGenerator::hyper_cube(triangulation, -3, 3);
    triangulation.refine_global(5);

    FE_Q<dim>       fe(1);
    DoFHandler<dim> dof_handler(triangulation);
    dof_handler.distribute_dofs(fe);

    Vector<double> values(dof_handler.n_dofs());

    VectorTools::interpolate(dof_handler, signed_distance, values);

    DataOut<dim> data_out;
    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(values, "signed_distance");

    data_out.build_patches();

    std::ofstream output(filename);
    data_out.write_vtu(output);
  }



  template <int dim>
  void
  test_infinite_cylinder_signed_distance()
  {
    deallog << "test_infinite_cylinder_signed_distance" << std::endl;

    const Point<dim> axis_point;

    Tensor<1, dim> axis;
    axis[dim - 1] = 1.0; // align with last coordinate direction

    const double radius = 1.0;

    const Functions::SignedDistance::InfiniteCylinder<dim> signed_distance(
      radius, axis, axis_point);


    output_signed_distance_vtu<dim>(signed_distance,
                                    "cylinder_dim" + std::to_string(dim) +
                                      ".vtu");
    Point<dim> p;

    // ---------------------------------------------------------
    deallog << "on axis (singular case)" << std::endl;
    p = Point<dim>();
    print_values_at_point(signed_distance, p);

    // ---------------------------------------------------------
    deallog << "inside" << std::endl;
    p    = Point<dim>();
    p[0] = 0.5; // r = 0.5 < 1
    print_values_at_point(signed_distance, p);

    // ---------------------------------------------------------
    deallog << "on surface" << std::endl;
    p    = Point<dim>();
    p[0] = 1.0; // r = 1
    print_values_at_point(signed_distance, p);

    // ---------------------------------------------------------
    deallog << "outside" << std::endl;
    p    = Point<dim>();
    p[0] = 2.0; // r = 2
    print_values_at_point(signed_distance, p);

    // ---------------------------------------------------------
    deallog << "diagonal in cross-section" << std::endl;
    p = Point<dim>();
    if (dim > 1)
      {
        p[0] = 1.0;
        p[1] = 1.0; // r = sqrt(2)
      }
    print_values_at_point(signed_distance, p);

    // ---------------------------------------------------------
    deallog << "shift along axis" << std::endl;
    p          = Point<dim>();
    p[0]       = 0.5;
    p[dim - 1] = 10.0; // move along axis
    print_values_at_point(signed_distance, p);

    // ---------------------------------------------------------
    deallog << "far outside" << std::endl;
    p    = Point<dim>();
    p[0] = 3.0;
    if (dim > 1)
      p[1] = 4.0; // r = 5
    print_values_at_point(signed_distance, p);

    deallog << "gradient tests" << std::endl;

    {
      p    = Point<dim>();
      p[0] = 1.0;

      deallog << "surface point in x-direction" << std::endl;
      print_gradient_at_point(signed_distance, p);
    }

    {
      p    = Point<dim>();
      p[0] = 2.0;

      deallog << "outside point in x-direction" << std::endl;
      print_gradient_at_point(signed_distance, p);
    }

    if (dim > 1)
      {
        p    = Point<dim>();
        p[0] = 1.0;
        p[1] = 1.0;

        deallog << "diagonal cross-section point" << std::endl;
        print_gradient_at_point(signed_distance, p);
      }

    {
      p          = Point<dim>();
      p[0]       = 1.0;
      p[dim - 1] = 10.0;

      deallog << "shifted along axis" << std::endl;
      print_gradient_at_point(signed_distance, p);
    }
  }



  template <int dim>
  void
  test_infinite_cylinder_signed_distance_shifted_axis()
  {
    deallog << "shifted + angled axis test" << std::endl;

    Point<dim> axis_point;
    axis_point[0] = 1.0;
    if (dim > 1)
      axis_point[1] = -1.0;

    Tensor<1, dim> axis;
    if (dim == 1)
      axis[0] = 1.0;
    else
      {
        axis[0] = 1.0;
        axis[1] = 1.0;
        if (dim > 2)
          axis[2] = 0.0;
        axis /= axis.norm();
      }

    const double radius = 1.0;

    const Functions::SignedDistance::InfiniteCylinder<dim> signed_distance(
      radius, axis, axis_point);

    output_signed_distance_vtu<dim>(signed_distance,
                                    "cylinder_shifted_dim" +
                                      std::to_string(dim) + ".vtu");

    // ---------------------------------------------------------
    deallog << "point on axis" << std::endl;
    Point<dim> p = axis_point;
    print_values_at_point(signed_distance, p);

    // ---------------------------------------------------------
    deallog << "point on surface (orthogonal direction)" << std::endl;

    Point<dim> n;
    if (dim > 1)
      {
        // construct a vector orthogonal to axis in 2D/3D
        n[0] = axis[1];
        n[1] = -axis[0];
        if (dim > 2)
          n[2] = 0.0;

        n /= n.norm();
      }
    else
      n[0] = 0.0;

    p = axis_point;
    for (unsigned int i = 0; i < dim; ++i)
      p[i] += radius * n[i];

    print_values_at_point(signed_distance, p);
    print_gradient_at_point(signed_distance, p);

    // ---------------------------------------------------------
    deallog << "inside point" << std::endl;

    p = axis_point;
    for (unsigned int i = 0; i < dim; ++i)
      p[i] += 0.5 * radius * n[i];

    print_values_at_point(signed_distance, p);

    // ---------------------------------------------------------
    deallog << "outside point" << std::endl;

    p = axis_point;
    for (unsigned int i = 0; i < dim; ++i)
      p[i] += 2.0 * radius * n[i];

    print_values_at_point(signed_distance, p);

    // ---------------------------------------------------------
    deallog << "shift along axis" << std::endl;

    p = axis_point;
    for (unsigned int i = 0; i < dim; ++i)
      p[i] += 0.5 * radius * n[i];

    for (unsigned int i = 0; i < dim; ++i)
      p[i] += 5.0 * axis[i]; // move along axis

    print_values_at_point(signed_distance, p);
    print_gradient_at_point(signed_distance, p);
  }



  template <int dim>
  void
  run_test()
  {
    deallog << "dim = " << dim << std::endl;
    deallog << std::endl;

    test_infinite_cylinder_signed_distance<dim>();
    test_infinite_cylinder_signed_distance_shifted_axis<dim>();

    deallog << std::endl;
  }
} // namespace



int
main()
{
  initlog();
  run_test<2>();
  run_test<3>();
}
