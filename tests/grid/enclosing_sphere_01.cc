// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2020 by the deal.II authors
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


// Computes a reasonably small enclosing ball on a variety of cells.
// The design of this test and part of the blessed file are taken
// from measure_et_al_02 test.


#include <deal.II/base/exceptions.h>
#include <deal.II/base/geometry_info.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <limits>

#include "../tests.h"

#define PRECISION 5


template <int dim>
void
create_triangulation(const unsigned int, Triangulation<dim> &)
{
  Assert(false, ExcNotImplemented());
}


template <>
void
create_triangulation(const unsigned int case_no, Triangulation<2> &tria)
{
  switch (case_no)
    {
      case 0:
        GridGenerator::hyper_cube(tria, 1., 3.);
        break;
      case 1:
        {
          GridGenerator::hyper_cube(tria, 1., 3.);
          Point<2> &v0 = tria.begin_active()->vertex(0);
          v0(0)        = 0.;
          break;
        }
      case 2:
        {
          GridGenerator::hyper_cube(tria, 1., 3.);
          Point<2> &v0 = tria.begin_active()->vertex(0);
          v0(0)        = 0.;
          Point<2> &v3 = tria.begin_active()->vertex(3);
          v3(0)        = 4.;
          break;
        }
      case 3:
        {
          GridGenerator::hyper_cube(tria, 1., 3.);
          Point<2> &v0 = tria.begin_active()->vertex(0);
          v0           = Point<2>(1.9, 1.9);

          Point<2> &v3 = tria.begin_active()->vertex(3);
          v3           = Point<2>(3.1, 3.);
          break;

          //
          //  y^  2-------3
          //   |  |       |
          //   |  |       |
          //   |  |       |
          //   |  0-------1
          //   *------------>x
          //
          //        |
          //        v
          //
          // vertices 0 and 3 are moved
          // such that the initial guess for the ball
          // (with its diameter as the largest diagonal (between vertices 1-2))
          // doesn't enclose vertex 3
          //
          //  y^  2--------+3
          //   |  |        |
          //   |  |        |
          //   |  | 0      |
          //   |  +--------1
          //   *------------>x
        }
      default:
        Assert(false, ExcNotImplemented());
    };
}


template <>
void
create_triangulation(const unsigned int case_no, Triangulation<3> &tria)
{
  switch (case_no)
    {
      case 0:
        GridGenerator::hyper_cube(tria, 1., 3.);
        break;
      case 1:
      case 2: // like case 1
        {
          GridGenerator::hyper_cube(tria, 1., 3.);
          Point<3> &v0 = tria.begin_active()->vertex(0);
          v0(0)        = 0.;
          break;
        }
      case 3:
        {
          //       6-------7        6-------7
          //      /|       |       /       /|
          //     / |       |      /       / |
          //    /  |       |     /       /  |
          //   4   |       |    4-------5   |
          //   |   2-------3    |       |   3
          //   |  /       /     |       |  /
          //   | /       /      |       | /
          //   |/       /       |       |/
          //   0-------1        0-------1
          //
          //
          //                |
          //                v
          //
          // moving vertices 1 and 6 such that diagonal 16 is smaller than all
          // other diagonals but vertex 6 lies outside all the balls constructed
          // using largest diagonals
          //
          //
          //      6                6
          //       +-------7        +-------7
          //      /|       |       /       /|
          //     / |       |      /       / |
          //    /  |       |     /       /  |
          //   4   |       |    4-------5   |
          //   |   +-------3    |       |   3
          //   |  /       /     |       |  /
          //   | /       /      |       | /
          //   |/     1 /       |       |/
          //   0-------+        0-------+
          //

          GridGenerator::hyper_cube(tria, 1., 3.);
          Point<3> &v1 = tria.begin_active()->vertex(1);
          Point<3> &v6 = tria.begin_active()->vertex(6);
          v1 += Point<3>(-0.9, 0.9, 0.9); // v1 was (3.,1.,1.)
          v6 += Point<3>(-0.2, 0.2, 0.2); // v7 was (1.,3.,3.)
          const Point<3> &v0 = tria.begin_active()->vertex(0);
          const Point<3> &v7 = tria.begin_active()->vertex(7);
          AssertThrow(v1.distance(v6) < v0.distance(v7),
                      ExcMessage("Vertices not moved as drawn above"));
          break;
        }
      default:
        Assert(false, ExcNotImplemented());
    };
}


template <int dim>
void
test()
{
  Triangulation<dim> tria;
  for (unsigned int case_no = 0; case_no < 4; ++case_no)
    {
      create_triangulation(case_no, tria);
      const std::pair<Point<dim>, double> smallest_sphere =
        tria.begin_active()->enclosing_ball();
      const double &    radius = smallest_sphere.second;
      const Point<dim> &center = smallest_sphere.first;

      deallog << "dim" << dim << ":case" << case_no << ":diameter=" << radius
              << ":center=" << center << std::endl;

      // Check that all the vertices are within the sphere
      // (sphere with thickness 100. *std::numeric_limits<double>::epsilon())
      for (const unsigned int v : GeometryInfo<dim>::vertex_indices())
        AssertThrow(std::fabs(center.distance(tria.begin_active()->vertex(v))) <
                      radius + 100. * std::numeric_limits<double>::epsilon(),
                    ExcInternalError());

      tria.clear();
    }
  deallog << "PASSED!" << std::endl;
}


int
main()
{
  initlog();
  deallog << std::setprecision(PRECISION);

  test<2>();
  test<3>();
}
