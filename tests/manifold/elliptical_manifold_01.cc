// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Generates a simple hyper_shell over an EllipticalManifold and
// test refinements.
// The triangulation produced in this test is part of the documentation
// of EllipticalManifold class as elliptical_hyper_shell.png.

#include "../tests.h"


// all include files you need here
#include <deal.II/base/exceptions.h>
#include <deal.II/base/numbers.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <vector>


// Anonymous namespace collecting a set of functions that creates an elliptical
// hyper_shell made of 8 elements.
namespace
{
  std::vector<Point<2>>
  generate_shell_points(const Point<2> &center,
                        const double    radius0,
                        const double    radius1,
                        const double    eccentricity)
  {
    const std::array<double, 2> radii{
      {radius0 * eccentricity, radius1 * eccentricity}};
    const std::array<double, 5> angle{
      {0, numbers::PI_4, numbers::PI_2, 3.0 * numbers::PI_4, numbers::PI}};
    std::vector<Point<2>> points;
    // Create an elliptical manifold to use push_forward()
    Tensor<1, 2> axis;
    axis[0] = 1.0;
    axis[1] = 0.0;
    EllipticalManifold<2, 2> mani(center, axis, eccentricity);
    for (auto &i : radii)
      {
        for (auto &j : angle)
          {
            points.emplace_back(mani.push_forward(Point<2>(i, j)) + center);
          }
      }
    std::array<int, 6> dup{{1, 2, 3, 6, 7, 8}};
    for (int &i : dup)
      {
        Point<2> pt = points[i];
        pt[1]       = -pt[1];
        points.emplace_back(pt);
      }
    return points;
  }



  // Generate an hyper_shell over an EllipticalManifold having an arbitrary
  // center, and the major axis oriented in the direction of the x-axis.
  //
  // inner_radius and outer_radius parameters correspond to the
  // semi-major axis length for the inner and outer ellipsis.
  //
  // Eccentricity must be in range ]0,1[
  //
  void
  build_simple_hyper_shell(Triangulation<2, 2> &grid,
                           const Point<2>      &center,
                           const double         inner_radius,
                           const double         outer_radius,
                           const double         eccentricity)
  {
    unsigned int          cell[][4] = {{5, 6, 0, 1},
                                       {6, 7, 1, 2},
                                       {8, 3, 7, 2},
                                       {9, 4, 8, 3},
                                       {5, 0, 13, 10},
                                       {13, 10, 14, 11},
                                       {15, 14, 12, 11},
                                       {9, 15, 4, 12}};
    std::vector<Point<2>> points =
      generate_shell_points(center, inner_radius, outer_radius, eccentricity);
    std::vector<CellData<2>> cells(8, CellData<2>());
    for (int i = 0; i < 8; ++i)
      for (int j = 0; j < 4; ++j)
        cells[i].vertices[j] = cell[i][j];
    SubCellData scdata;
    grid.create_triangulation(points, cells, scdata);
    // assign manifold
    Tensor<1, 2> axis;
    axis[0] = 1.0;
    axis[1] = 0.0;
    grid.set_manifold(0, EllipticalManifold<2, 2>(center, axis, eccentricity));
    grid.set_all_manifold_ids(0);
  }



} // namespace



// Helper function
// Generate a simple hyper_shell over an elliptical manifold centered at the
// origin. Major axis is the x-axis.
template <int dim, int spacedim>
void
test(unsigned int ref = 1)
{
  Triangulation<dim, spacedim> tria;
  switch (dim)
    {
      case 2:
        build_simple_hyper_shell(tria, Point<spacedim>(), 1.0, 2.0, 0.8);
        break;
      default:
        AssertThrow(0, ExcNotImplemented());
        break;
    }

  tria.refine_global(3);

  GridOut gridout;
  gridout.write_msh(tria, deallog.get_file_stream());
}



int
main()
{
  initlog();
  test<2, 2>();
}
