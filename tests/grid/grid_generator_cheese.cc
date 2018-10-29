// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2017 by the deal.II authors
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

// Test output for GridGenerator::cheese()

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

  std::vector<unsigned int> holes(d);
  holes[0] = 3;
  holes[1] = 2;
  GridGenerator::cheese(tr, holes);

  GridOut gout;
  gout.write_vtk(tr, os);
}

void
dim_3(std::ostream &os)
{
  const unsigned int d = 3;
  Triangulation<d>   tr;

  std::vector<unsigned int> holes(d);
  holes[0] = 3;
  holes[1] = 2;
  holes[2] = 4;
  GridGenerator::cheese(tr, holes);

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
