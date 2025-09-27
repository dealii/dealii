// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test output for GridGenerator::simplex()

#include <deal.II/base/tensor.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"



void
dim_2(std::ostream &os)
{
  const unsigned int d = 2;
  Triangulation<d>   tr;

  std::vector<Point<d>> vertices(d + 1);
  vertices[1][0] = 0.5;
  vertices[1][1] = .85;
  vertices[2][0] = -0.5;
  vertices[2][1] = .85;
  GridGenerator::simplex(tr, vertices);

  GridOut gout;
  gout.write_vtk(tr, os);
}

void
dim_3(std::ostream &os)
{
  const unsigned int d = 3;
  Triangulation<d>   tr;

  std::vector<Point<d>> vertices(d + 1);
  vertices[0][0] = 1.;
  vertices[0][1] = 0.;
  vertices[0][2] = .7;
  vertices[1][0] = -1.;
  vertices[1][1] = 0.;
  vertices[1][2] = .7;
  vertices[2][0] = 0.;
  vertices[2][1] = 1.;
  vertices[2][2] = -.7;
  vertices[3][0] = 0.;
  vertices[3][1] = -1.;
  vertices[3][2] = -.7;
  GridGenerator::simplex(tr, vertices);

  GridOut gout;
  gout.write_vtk(tr, os);
}


int
main()
{
  initlog(true);
  std::ostream &logfile = deallog.get_file_stream();
  dim_2(logfile);
  dim_3(logfile);
}
