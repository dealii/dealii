// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2012 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// take a disconnected 2d mesh and check that we can find an arbitrary point's
// cell in it. We consider a special triangulation, where the point p does not
// lie in a cell adjacent to the vertex with minimal distance to p. The test
// should fail for all revisions <= 25704M.

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"


void
create_coarse_grid(Triangulation<2> &coarse_grid)
{
  static const Point<2> vertices_1[] = {
    Point<2>(0., 0.), // 0
    Point<2>(1., 0.), // 1
    Point<2>(1., 1.), // 2
    Point<2>(0., 1.), // 3

    Point<2>(1.1, 0.),      // 4
    Point<2>(1.1, 1. / 2.), // 5
    Point<2>(1.3, 0.),      // 6
    Point<2>(1.3, 1. / 2.), // 7
  };
  const unsigned int n_vertices = sizeof(vertices_1) / sizeof(vertices_1[0]);

  const std::vector<Point<2>> vertices(&vertices_1[0], &vertices_1[n_vertices]);

  static const int cell_vertices[][GeometryInfo<2>::vertices_per_cell] = {
    {0, 1, 3, 2},
    {4, 6, 5, 7},
  };
  const unsigned int n_cells = sizeof(cell_vertices) / sizeof(cell_vertices[0]);

  std::vector<CellData<2>> cells(n_cells, CellData<2>());
  for (unsigned int i = 0; i < n_cells; ++i)
    {
      for (const unsigned int j : GeometryInfo<2>::vertex_indices())
        cells[i].vertices[j] = cell_vertices[i][j];
      cells[i].material_id = 0;
    }

  coarse_grid.create_triangulation(vertices, cells, SubCellData());
}


void
check(Triangulation<2> &tria)
{
  Point<2> p(0.99, 1. / 2.);

  Triangulation<2>::active_cell_iterator cell =
    GridTools::find_active_cell_around_point(tria, p);

  deallog << cell << std::endl;
  for (const unsigned int v : GeometryInfo<2>::vertex_indices())
    deallog << '<' << cell->vertex(v) << "> ";
  deallog << std::endl;

  AssertThrow(p.distance(cell->center()) < cell->diameter() / 2,
              ExcInternalError());
}


int
main()
{
  initlog();

  {
    Triangulation<2> coarse_grid;
    create_coarse_grid(coarse_grid);
    check(coarse_grid);
  }
}
