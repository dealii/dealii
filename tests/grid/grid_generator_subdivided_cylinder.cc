// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2016 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------



// Test GridTools::consistently_order_cells

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
