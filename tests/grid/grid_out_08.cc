// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test GridOut::write_vtu for 1d meshes in 3d

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

  GridOut           grid_out;
  GridOutFlags::Vtu vtu_flags;
  vtu_flags.compression_level = DataOutBase::CompressionLevel::best_compression;
  grid_out.set_flags(vtu_flags);
  grid_out.write_vtu(tria, deallog.get_file_stream());
}


int
main()
{
  initlog();
  deallog.get_file_stream() << std::setprecision(2);

  test();
}
