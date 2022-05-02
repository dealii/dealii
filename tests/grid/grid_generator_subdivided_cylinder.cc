// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2021 by the deal.II authors
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



// Test GridReordering<3>::reorder_cells

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"



void
dim_3(std::ostream &os)
{
  const unsigned int d = 3;
  Triangulation<d>   tr;

  GridGenerator::subdivided_cylinder(tr, 4, 1, 1);

  GridOut gout;
  gout.write_vtk(tr, os);
}

int
main()
{
  initlog(true);
  std::ostream &logfile = deallog.get_file_stream();
  dim_3(logfile);
}
