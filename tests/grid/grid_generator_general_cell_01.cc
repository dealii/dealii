// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test output for GridGenerator::general_cell()

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"

void
dim_2(std::ostream &os)
{
  std::vector<Point<2>> vertices(4);
  vertices[0][0] = -1.;
  vertices[0][1] = -1.;
  vertices[1][0] = 1.;
  vertices[1][1] = -1.5;
  vertices[2][0] = 1.5;
  vertices[2][1] = 1.5;
  vertices[3][0] = 2.;
  vertices[3][1] = 0.5;

  Triangulation<2> tria;
  GridGenerator::general_cell<2>(tria, vertices);

  GridOut gout;
  gout.write_vtk(tria, os);
}

void
dim_3(std::ostream &os)
{
  std::vector<Point<3>> vertices(8);
  vertices[0][0] = -1.;
  vertices[0][1] = -1.;
  vertices[0][2] = -1.;
  vertices[1][0] = 1.;
  vertices[1][1] = -1.5;
  vertices[1][2] = -1.5;
  vertices[2][0] = 2.;
  vertices[2][1] = 1.5;
  vertices[2][2] = -2.;
  vertices[3][0] = 2.5;
  vertices[3][1] = 0.5;
  vertices[3][2] = -3.;
  vertices[4][0] = -1.;
  vertices[4][1] = -1.;
  vertices[4][2] = 1.;
  vertices[5][0] = 1.;
  vertices[5][1] = -1.5;
  vertices[5][2] = 1.5;
  vertices[6][0] = 2.;
  vertices[6][1] = 1.5;
  vertices[6][2] = 2.;
  vertices[7][0] = 2.;
  vertices[7][1] = 0.5;
  vertices[7][2] = 3.;

  Triangulation<3> tria;
  GridGenerator::general_cell<3>(tria, vertices);

  GridOut gout;
  gout.write_vtk(tria, os);
}


int
main()
{
  initlog(true);
  std::ostream &logfile = deallog.get_file_stream();
  dim_2(logfile);
  dim_3(logfile);
}
