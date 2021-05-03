// ---------------------------------------------------------------------
//
// Copyright (C) 2007 - 2018 by the deal.II authors
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

/**
 * Test the classes Functions::LevelSet::Sphere and Functions::LevelSet::Plane,
 * by creating an object of each class and calling the value, gradient and
 * hessian functions at a point where we know what the function values should
 * be.
 */

#include <deal.II/base/function_level_set.h>

#include "../tests.h"


namespace
{
  template <int dim>
  void
  print_derivatives_at_point(const Function<dim> &function,
                             const Point<dim> &   point)
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

    const Functions::LevelSet::Sphere<dim> level_set(center, radius);

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

    const Functions::LevelSet::Plane<dim> level_set(point_in_plane, normal);

    Point<dim> point;
    for (unsigned int i = 0; i < dim; i++)
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
