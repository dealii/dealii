// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 by the deal.II authors
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
 * Test the classes Functions::SignedDistance::ZalesakDisk by creating an
 * object of each class and calling the value functions at a point where we know
 * what the function values should be.
 */

#include <deal.II/base/function_signed_distance.h>

#include "../tests.h"


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
  test_disk_level_set()
  {
    deallog << "test_zalesak_disk_level_set" << std::endl;

    const Point<dim> center;
    const double     radius = 1;

    const Functions::SignedDistance::ZalesakDisk<dim> level_set(center,
                                                                radius,
                                                                radius * 0.5,
                                                                0.75 * radius);

    {
      Point<dim> p;

      deallog << "center" << std::endl;
      print_values_at_point(level_set, p);
      deallog << "inside" << std::endl;
      p[0] += 0.1;
      print_values_at_point(level_set, p);
      p[0] = 0.3;
      print_values_at_point(level_set, p);
    }
    deallog << "notch" << std::endl;
    Point<dim> p;
    p[dim - 1]              = -0.45;
    p[std::max(0, dim - 2)] = -0.1;
    print_values_at_point(level_set, p);
    deallog << "outside" << std::endl;
    p[dim - 1] = -2;
    print_values_at_point(level_set, p);
    p[dim - 1]              = -1.2;
    p[0]                    = -0.3;
    p[std::max(0, dim - 2)] = -0.3;
    print_values_at_point(level_set, p);
    p[dim - 1]              = -2;
    p[0]                    = 0;
    p[std::max(0, dim - 2)] = 0;
    print_values_at_point(level_set, p);
    p[dim - 1] = 2;
    print_values_at_point(level_set, p);
    p[0] = -2;
    print_values_at_point(level_set, p);
    p[0] = 2;
    print_values_at_point(level_set, p);
  }

  template <int dim>
  void
  run_test()
  {
    deallog << "dim = " << dim << std::endl;
    deallog << std::endl;
    test_disk_level_set<dim>();
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
