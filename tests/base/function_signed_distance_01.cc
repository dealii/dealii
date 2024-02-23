// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

/**
 * Test the classes Functions::SignedDistance::Sphere and
 * Functions::SignedDistance::Plane, by creating an object of each class and
 * calling the value, gradient and hessian functions at a point where we know
 * what the function values should be.
 */

#include <deal.II/base/function_signed_distance.h>

#include "../tests.h"


namespace
{
  template <int dim>
  void
  print_derivatives_at_point(const Function<dim> &function,
                             const Point<dim>    &point)
  {
    deallog << "point = " << point << std::endl;
    deallog << "value = " << function.value(point) << std::endl;
    deallog << "gradient = " << function.gradient(point) << std::endl;
    deallog << "Hessian = " << function.hessian(point) << std::endl;
  }



  template <int dim>
  void
  test_sphere_level_set()
  {
    deallog << "test_sphere_level_set" << std::endl;

    const Point<dim> center;
    const double     radius = 1;

    const Functions::SignedDistance::Sphere<dim> level_set(center, radius);

    Point<dim> point;
    point[0] = 2 * radius;

    print_derivatives_at_point(level_set, point);
  }



  template <int dim>
  void
  test_plane_level_set()
  {
    deallog << "test_plane_level_set" << std::endl;

    const Point<dim>     point_in_plane;
    const Tensor<1, dim> normal = Point<dim>::unit_vector(0);

    const Functions::SignedDistance::Plane<dim> level_set(point_in_plane,
                                                          normal);

    Point<dim> point;
    for (unsigned int i = 0; i < dim; ++i)
      point[i] = 1;

    print_derivatives_at_point(level_set, point);
  }



  template <int dim>
  void
  run_test()
  {
    deallog << "dim = " << dim << std::endl;
    deallog << std::endl;
    test_sphere_level_set<dim>();
    deallog << std::endl;
    test_plane_level_set<dim>();
    deallog << std::endl;
  }
} // namespace



int
main()
{
  initlog();
  run_test<1>();
  run_test<2>();
  run_test<3>();
}
