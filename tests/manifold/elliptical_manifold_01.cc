// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2018 by the deal.II authors
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
  /* The following system of equations
  // /
  // |x = ellPar*cosh(u)*cos(v)
  // |y = ellPar*sinh(u)*sin(v)
  // \
  // allows to transform a point (u,v) in the chart space to a
  // point (x,y) in the cartesian space.
  //
  // Setting v=0, the system of equations becomes
  // /
  // |x = x_0 = ellPar*cosh(u)
  // |y = 0
  // \
  // which represents the crossing point of an ellipsis, non-rotated and
  // centered at the origin of the cartesian system, with its major axis.
  //
  // For the sake of completeness, by setting v=pi/2, one finds y_0 =
  // ellPar*sinh(u) as the crossing point on the minor axis.
  //
  // This function takes in input a point pt={x_0,v}
  // and transforms it into (x,y) in accordance with the first system of
  // equations above.
  //
  // On a final note, the condition x_0 >= ellPar has to be satisfied. This is
  // required because min(cosh(u)) = 1 is obtained for u = 0, and imposing v = 0
  // the minimum value that x can assume is x = ellPar.
  // For the sake of completeness, it's worth to notice that there would be no
  // restrictions on the input if we passed y_0, instead of x_0, as pt[0].
  */
  Point<2>
  chart_to_cartesian(const Point<2> &pt, const double ellPar)
  {
    Point<2>     c;
    const double ch = pt[0] / ellPar;
    const double sh = std::sqrt(ch * ch - 1.0);
    c[0]            = ellPar * ch * std::cos(pt[1]);
    c[1]            = ellPar * sh * std::sin(pt[1]);
    return c;
  }



  std::vector<Point<2>>
  generate_shell_points(const Point<2> &center,
                        const double    radius0,
                        const double    radius1,
                        const double    ellPar)
  {
    const std::array<double, 2> radii{{radius0, radius1}};
    const std::array<double, 5> angle{
      {0, numbers::PI_4, numbers::PI_2, 3.0 * numbers::PI_4, numbers::PI}};
    std::vector<Point<2>> points;
    for (auto &i : radii)
      {
        for (auto &j : angle)
          {
            points.emplace_back(chart_to_cartesian(Point<2>(i, j), ellPar) +
                                center);
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



  // Generate an hyper_shell over an EllipticalManifold having an abitrary
  // center, and the major axis oriented in the direction of the x-axis.
  //
  // inner_radius and outer_radius parameters correspond to the
  // distances from the center of the manifold of the
  // crossing points of two concentric ellipsis with their major axis.
  //
  // Input parameters must respect the following constrains:
  // c_parameter < inner_radius < outer_radius.
  //
  // To understand the constrain on the c_parameter refer to the
  // documentation of function chart_to_cartesian() whithin this namespace.
  void build_simple_hyper_shell(Triangulation<2, 2> &grid,
                                const Point<2> &     center,
                                const double         inner_radius,
                                const double         outer_radius,
                                const double         c_param)
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
      generate_shell_points(center, inner_radius, outer_radius, c_param);
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
    grid.set_manifold(0, EllipticalManifold<2, 2>(center, axis, c_param));
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
  // deallog << "Testing dim " << dim << ", spacedim " << spacedim << std::endl;

  Triangulation<dim, spacedim> tria;
  switch (dim)
    {
      case 2:
        build_simple_hyper_shell(tria, Point<spacedim>(), 1.5, 1.8, 1.4);
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
