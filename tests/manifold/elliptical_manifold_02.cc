// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2019 by the deal.II authors
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

// This test verifies the correct functioning of pull_back() and push_forward().
// It checks that pull_back(push_forward(input)) returns input within a
// given tolerance.

#include "../tests.h"


// all include files you need here
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#define TOLERANCE 1e-14



// Helper functions
template <int dim, int spacedim>
Point<spacedim>
test_push(const Tensor<1, spacedim> &axis,
          const double               eccentricity,
          const Point<spacedim> &    space_point,
          unsigned int               ref = 1)
{
  EllipticalManifold<dim, spacedim> manifold(Point<spacedim>(),
                                             axis,
                                             eccentricity);
  const Point<spacedim> chart_point(manifold.push_forward(space_point));
  deallog << space_point << " -> ";
  deallog << chart_point << std::endl;
  return chart_point;
}



template <int dim, int spacedim>
Point<spacedim>
test_pull(const Tensor<1, spacedim> &axis,
          const double               eccentricity,
          const Point<spacedim> &    chart_point,
          unsigned int               ref = 1)
{
  EllipticalManifold<dim, spacedim> manifold(Point<spacedim>(),
                                             axis,
                                             eccentricity);
  const Point<spacedim> space_point(manifold.pull_back(chart_point));
  deallog << space_point << " <- ";
  deallog << chart_point << std::endl;
  return space_point;
}



// Function that tests pull_back() and push_forward().
void
local_test(const Tensor<1, 2> &axis,
           const double        eccentricity,
           const Point<2> &    space_point)
{
  const Point<2> pt1 = test_push<2, 2>(axis, eccentricity, space_point);
  const Point<2> pt2 = test_pull<2, 2>(axis, eccentricity, pt1);
  if ((space_point - pt2).norm() < TOLERANCE)
    {
      deallog << "OK" << std::endl;
    }
  else
    {
      deallog << "FAILED" << std::endl;
    }
}



void
test()
{
  {
    // test on a default manifold
    const Tensor<1, 2> axis({1.0, 0.0});
    const double       eccentricity(0.1);
    local_test(axis, eccentricity, Point<2>(1.1, 0.0));
    local_test(axis, eccentricity, Point<2>(2.0, 0.0));
    local_test(axis, eccentricity, Point<2>(3.0, 0.0));
    local_test(axis, eccentricity, Point<2>(4.0, 2.0));
  }
  {
    // test on a rotated manifold
    const Tensor<1, 2> axis({1.0, -1.0});
    const double       eccentricity(0.1);
    local_test(axis, eccentricity, Point<2>(1.1, 0.25 * TOLERANCE));
    local_test(axis, eccentricity, Point<2>(2.0, 0.25 * TOLERANCE));
    local_test(axis, eccentricity, Point<2>(3.0, 0.25 * TOLERANCE));
    local_test(axis, eccentricity, Point<2>(4.0, 2.0));
  }
}



int
main()
{
  initlog();
  test();
}
