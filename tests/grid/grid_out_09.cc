// ---------------------------------------------------------------------
//
// Copyright (C) 2002 - 2020 by the deal.II authors
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


// test GridOut::write_gnuplot for 1d meshes in 3d

#include <deal.II/base/geometry_info.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"

void
test()
{
  Triangulation<1, 3> tria;

  const unsigned int    dim          = 1;
  static const Point<3> vertices_1[] = {
    Point<3>(0, 0, 0), Point<3>(1, 0, 0), Point<3>(1, 0, 0), Point<3>(1, 1, 0),
    Point<3>(1, 1, 0), Point<3>(0, 1, 0), Point<3>(0, 1, 0), Point<3>(0, 0, 0),
    Point<3>(0, 0, 1), Point<3>(1, 0, 1), Point<3>(1, 0, 1), Point<3>(1, 1, 1),
    Point<3>(1, 1, 1), Point<3>(0, 1, 1), Point<3>(0, 1, 1), Point<3>(0, 0, 1),
    Point<3>(0, 0, 0), Point<3>(0, 0, 1), Point<3>(1, 0, 0), Point<3>(1, 0, 1),
    Point<3>(1, 1, 0), Point<3>(1, 1, 1), Point<3>(0, 1, 0), Point<3>(0, 1, 1)};
  const unsigned int n_vertices = sizeof(vertices_1) / sizeof(vertices_1[0]);
  const std::vector<Point<3>> vertices(&vertices_1[0], &vertices_1[n_vertices]);

  static const int cell_vertices[][GeometryInfo<dim>::vertices_per_cell] = {
    {0, 1},
    {2, 3},
    {4, 5},
    {6, 7},
    {8, 9},
    {10, 11},
    {12, 13},
    {14, 15},
    {16, 17},
    {18, 19},
    {20, 21},
    {22, 23}};
  const unsigned int n_cells = sizeof(cell_vertices) / sizeof(cell_vertices[0]);
  std::vector<CellData<dim>> cells(n_cells, CellData<dim>());

  for (unsigned int i = 0; i < n_cells; ++i)
    {
      for (const unsigned int j : GeometryInfo<dim>::vertex_indices())
        cells[i].vertices[j] = cell_vertices[i][j];
      cells[i].material_id = 0;
    }
  tria.create_triangulation(vertices, cells, SubCellData());

  tria.refine_global(1);

  GridOut grid_out;
  grid_out.write_gnuplot(tria, deallog.get_file_stream());
}


int
main()
{
  initlog();
  deallog << std::setprecision(2);

  test();
}
