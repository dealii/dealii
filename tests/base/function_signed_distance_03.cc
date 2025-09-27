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
 * Test the class Functions::SignedDistance::Rectangle by creating an object of
 * the class and calling the value function at points where we know what the
 * function values should be.
 */

#include <deal.II/base/function_signed_distance.h>

#include "../tests.h"

namespace
{
  template <int dim>
  void
  print_value_at_point(const Function<dim> &function, const Point<dim> &point)
  {
    deallog << "point = " << point << std::endl;
    deallog << "value = " << function.value(point) << std::endl;
  }


  void
  test_rectangle_signed_distance_1d()
  {
    static constexpr int dim = 1;

    deallog << "test_rectangle_signed_distance" << std::endl;
    deallog << "dim = " << dim << std::endl;
    deallog << std::endl;

    const Point<dim> bl(0);
    const Point<dim> tr(1);
    const Point<dim> center(0.5);

    const Functions::SignedDistance::Rectangle<dim> rectangle(bl, tr);

    deallog << "center" << std::endl;
    print_value_at_point(rectangle, center);
    deallog << "inside" << std::endl;
    print_value_at_point(rectangle, Point<dim>(0.25));
    print_value_at_point(rectangle, Point<dim>(0.75));
    deallog << "on rectangle" << std::endl;
    print_value_at_point(rectangle, bl);
    print_value_at_point(rectangle, tr);
    deallog << "outside" << std::endl;
    print_value_at_point(rectangle, Point<dim>(-0.5));
    print_value_at_point(rectangle, Point<dim>(1.5));

    deallog << std::endl;
  }



  void
  test_rectangle_signed_distance_2d()
  {
    static constexpr int dim = 2;

    deallog << "test_rectangle_signed_distance" << std::endl;
    deallog << "dim = " << dim << std::endl;
    deallog << std::endl;

    const Point<dim> bl(0, 0);
    const Point<dim> tr(1, 1);
    const Point<dim> center(0.5, 0.5);

    const Functions::SignedDistance::Rectangle<dim> rectangle(bl, tr);

    deallog << "center" << std::endl;
    print_value_at_point(rectangle, center);
    deallog << "inside" << std::endl;
    print_value_at_point(rectangle, Point<dim>(0.25, 0.25));
    print_value_at_point(rectangle, Point<dim>(0.75, 0.25));
    print_value_at_point(rectangle, Point<dim>(0.75, 0.75));
    print_value_at_point(rectangle, Point<dim>(0.25, 0.75));
    deallog << "on rectangle" << std::endl;
    print_value_at_point(rectangle, Point<dim>(0., 0.5));
    print_value_at_point(rectangle, Point<dim>(1.0, 0.5));
    print_value_at_point(rectangle, Point<dim>(1.0, 0.5123));
    print_value_at_point(rectangle, Point<dim>(0.0, 0.5123));
    deallog << "outside" << std::endl;

    for (double step = 0; step < 2.1; step += 0.5)
      {
        print_value_at_point(rectangle, Point<dim>(step, -1));
        print_value_at_point(rectangle, Point<dim>(step, 2));
      }
    for (double step = 0; step < 2.1; step += 0.5)
      {
        print_value_at_point(rectangle, Point<dim>(-1, step));
        print_value_at_point(rectangle, Point<dim>(2, step));
      }

    deallog << std::endl;
  }

  void
  test_rectangle_signed_distance_3d()
  {
    static constexpr int dim = 3;

    deallog << "test_rectangle_signed_distance" << std::endl;
    deallog << "dim = " << dim << std::endl;
    deallog << std::endl;

    const Point<dim> bl(0, 0, 0);
    const Point<dim> tr(1, 1, 1);
    const Point<dim> center(0.5, 0.5, 0.5);

    const Functions::SignedDistance::Rectangle<dim> rectangle(bl, tr);

    deallog << "center" << std::endl;
    print_value_at_point(rectangle, center);
    deallog << "inside" << std::endl;
    print_value_at_point(rectangle, Point<dim>(0.25, 0.25, 0.25));
    print_value_at_point(rectangle, Point<dim>(0.75, 0.25, 0.25));
    print_value_at_point(rectangle, Point<dim>(0.75, 0.75, 0.25));
    print_value_at_point(rectangle, Point<dim>(0.25, 0.75, 0.25));
    deallog << "on rectangle" << std::endl;
    print_value_at_point(rectangle, Point<dim>(0., 0.5, 0.5));
    print_value_at_point(rectangle, Point<dim>(0.5, 0., 0.5));
    print_value_at_point(rectangle, Point<dim>(0.5, 0.5, 0.));
    print_value_at_point(rectangle, Point<dim>(1., 0.5, 0.5));
    print_value_at_point(rectangle, Point<dim>(0.5, 1., 0.5));
    print_value_at_point(rectangle, Point<dim>(0.5, 0.5, 1.));
    deallog << "outside" << std::endl;
    print_value_at_point(rectangle, Point<dim>(-1, 0, 0));
    print_value_at_point(rectangle, Point<dim>(2, 0, 0));
    print_value_at_point(rectangle, Point<dim>(0, -1, 0));
    print_value_at_point(rectangle, Point<dim>(0, 2, 0));
    print_value_at_point(rectangle, Point<dim>(0, 0, -1));
    print_value_at_point(rectangle, Point<dim>(0, 0, 2));
    print_value_at_point(rectangle, Point<dim>(0, -1, 2));
    print_value_at_point(rectangle, Point<dim>(-1, 2, 0));
    print_value_at_point(rectangle, Point<dim>(2, 0, -1));
    deallog << std::endl;
  }
} // namespace



int
main()
{
  initlog();
  test_rectangle_signed_distance_1d();
  test_rectangle_signed_distance_2d();
  test_rectangle_signed_distance_3d();
}
