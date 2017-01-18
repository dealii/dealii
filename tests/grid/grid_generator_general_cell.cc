// ---------------------------------------------------------------------
//
// Copyright (C) 2016 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

// Test output for GridGenerator::general_cell()

#include "../tests.h"
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_out.h>

void dim2(std::ostream &os)
{
  std::vector<Point<2> > vertices(4);
  vertices[0](0) = -1.;
  vertices[0](1) = -1.;
  vertices[1](0) = 1.;
  vertices[1](1) = -1.5;
  vertices[2](0) = 1.5;
  vertices[2](1) = 1.5;
  vertices[3](0) = 2.;
  vertices[3](1) = 0.5;

  Triangulation<2> tria;
  GridGenerator::general_cell<2> (tria, vertices);

  GridOut gout;
  gout.write_vtk(tria, os);
}

void dim3(std::ostream &os)
{
  std::vector<Point<3> > vertices(8);
  vertices[0](0) = -1.;
  vertices[0](1) = -1.;
  vertices[0](2) = -1.;
  vertices[1](0) = 1.;
  vertices[1](1) = -1.5;
  vertices[1](2) = -1.5;
  vertices[2](0) = 2.;
  vertices[2](1) = 1.5;
  vertices[2](2) = -2.;
  vertices[3](0) = 2.5;
  vertices[3](1) = 0.5;
  vertices[3](2) = -3.;
  vertices[4](0) = -1.;
  vertices[4](1) = -1.;
  vertices[4](2) = 1.;
  vertices[5](0) = 1.;
  vertices[5](1) = -1.5;
  vertices[5](2) = 1.5;
  vertices[6](0) = 2.;
  vertices[6](1) = 1.5;
  vertices[6](2) = 2.;
  vertices[7](0) = 2.;
  vertices[7](1) = 0.5;
  vertices[7](2) = 3.;

  Triangulation<3> tria;
  GridGenerator::general_cell<3> (tria, vertices);

  GridOut gout;
  gout.write_vtk(tria, os);
}


int main()
{
  initlog(true);
  std::ostream &logfile = deallog.get_file_stream();
  dim2(logfile);
  dim3(logfile);
}
