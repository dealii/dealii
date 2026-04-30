// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2013 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

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

  // The mesh has 3*2=6 holes among (3*2+1)*(2*2+1)=7*5=35 possible
  // cell locations. So there must be 35-6=29 cells:
  Assert(tr.n_active_cells() == 29, ExcInternalError());

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

  // The mesh has 3*2*4=24 holes among
  // (3*2+1)*(2*2+1)*(4*2+1)=7*5*9=315 possible cell locations. So
  // there must be 315-24=291 cells:
  Assert(tr.n_active_cells() == 291, ExcInternalError());

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
