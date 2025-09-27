// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2006 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test GridTools::invert_cells_with_negative_measure

/*
  vertex layout:

  1   3   5

  0   2   4
 */

#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>

#include <string>

#include "../tests.h"

void
test()
{
  const static int        dim = 2;
  std::vector<Point<dim>> vertices(6);
  vertices[1][1] = 1;
  vertices[2][0] = 1;
  vertices[3][0] = 1;
  vertices[3][1] = 1;
  vertices[4][0] = 2;
  vertices[5][0] = 2;
  vertices[5][1] = 1;

  std::vector<CellData<dim>> cells(2);
  for (const unsigned int i : GeometryInfo<dim>::vertex_indices())
    cells[0].vertices[i] = i;
  cells[1].vertices[0] = 2;
  cells[1].vertices[1] = 4;
  cells[1].vertices[2] = 3;
  cells[1].vertices[3] = 5;

  SubCellData       subcelldata;
  const std::size_t n_cells_inverted =
    GridTools::invert_cells_with_negative_measure<>(vertices, cells);

  deallog << "We inverted " << n_cells_inverted << " cells." << std::endl;

  std::swap(cells[0].vertices[2], cells[0].vertices[3]);
  std::swap(cells[1].vertices[2], cells[1].vertices[3]);

  Triangulation<dim> tria;
  tria.create_triangulation(vertices, cells, subcelldata);

  std::ostream &logfile = deallog.get_file_stream();
  logfile << "---------------------------------------------" << std::endl
          << std::endl
          << std::endl;

  GridOut grid_out;
  grid_out.set_flags(GridOutFlags::Ucd(true));
  grid_out.write_ucd(tria, logfile);
}

int
main()
{
  initlog(false, std::ios_base::fmtflags());
  test();
}
