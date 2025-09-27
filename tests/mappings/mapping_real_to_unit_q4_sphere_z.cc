// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2012 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// something like the original of mapping_real_to_unit_q4_sphere: we
// wanted to know whether a point is inside or outside, but got an
// exception before r25651

#include <deal.II/base/utilities.h>

#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"


void
test_real_to_unit_cell()
{
  const unsigned int dim = 3;

  // define a spherical cap boundary
  // to be used for one of the faces
  const double radius =
    Point<dim>(1.43757e-10, 4.48023e+06, 4.48023e+06).norm();
  SphericalManifold<dim> boundary;

  // create the mesh: a single cell
  // with the following coordinates:
  std::vector<Point<dim>> vertices;
  vertices.push_back(Point<dim>(6.70384e-11, 3.17728e+06, 3.17728e+06));
  vertices.push_back(Point<dim>(-1.46060e+06, 3.99043e+06, 1.46060e+06));
  vertices.push_back(Point<dim>(-1.46060e+06, 1.46060e+06, 3.99043e+06));
  vertices.push_back(Point<dim>(-2.59424e+06, 2.59424e+06, 2.59424e+06));
  vertices.push_back(Point<dim>(1.43757e-10, 4.48023e+06, 4.48023e+06));
  vertices.push_back(Point<dim>(-2.05956e+06, 5.62684e+06, 2.05956e+06));
  vertices.push_back(Point<dim>(-2.05956e+06, 2.05956e+06, 5.62684e+06));
  vertices.push_back(Point<dim>(-3.65809e+06, 3.65809e+06, 3.65809e+06));
  // the points above don't show
  // enough digits to have the same
  // outer radius, so normalize the
  // four outer ones
  for (unsigned int v = 4; v < 8; ++v)
    vertices[v] *= radius / vertices[v].norm();
  std::vector<CellData<dim>> cells;
  {
    CellData<dim> d;
    for (unsigned int i = 0; i < 8; ++i)
      d.vertices[i] = i;
    cells.push_back(d);
  }
  Triangulation<dim> triangulation;
  triangulation.create_triangulation(vertices, cells, SubCellData());

  // set the boundary indicator for
  // one face and adjacent edges of
  // the single cell
  triangulation.set_manifold(1, boundary);
  triangulation.begin_active()->face(5)->set_all_boundary_ids(1);

  // now see if the point is inside
  // or outside
  const Point<dim> p(-3.56413e+06, 1.74215e+06, 2.14762e+06);
  Assert(triangulation.begin_active()->point_inside(p) == false,
         ExcInternalError());
  deallog << "Point is outside!" << std::endl;
}


int
main()
{
  initlog();

  test_real_to_unit_cell();
}
