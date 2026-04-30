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

// Test the variant of GridGenerator::cheese() that takes a voxel
// mask. For simplicity, use the same set-up as for the
// grid_generator_cheese.cc test where the voxels are evenly spaced in
// the dim_2() and dim_3() functions. In dim_2a(), use a 4x4 mesh in
// which the middle two 2x2 are empty; this create an unused vertex
// that the function internally has to remove.

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

  const bool pixel_mask[7][5] = {{true, true, true, true, true},
                                 {true, false, true, false, true},
                                 {true, true, true, true, true},
                                 {true, false, true, false, true},
                                 {true, true, true, true, true},
                                 {true, false, true, false, true},
                                 {true, true, true, true, true}};
  GridGenerator::cheese(tr, Table<2, bool>(7, 5, &pixel_mask[0][0]), 1.);

  // The mesh has 3*2=6 holes among (3*2+1)*(2*2+1)=7*5=35 possible
  // cell locations. So there must be 35-6=29 cells:
  AssertDimension(tr.n_active_cells(), 29);

  GridOut gout;
  gout.write_vtk(tr, os);
}



void
dim_2a(std::ostream &os)
{
  const unsigned int d = 2;
  Triangulation<d>   tr;

  const bool pixel_mask[4][4] = {{true, true, true, true},
                                 {true, false, false, true},
                                 {true, false, false, true},
                                 {true, true, true, true}};
  GridGenerator::cheese(tr, Table<2, bool>(4, 4, &pixel_mask[0][0]), 2);

  AssertDimension(tr.n_active_cells(), 12);

  GridOut gout;
  gout.write_vtk(tr, os);
}



void
dim_3(std::ostream &os)
{
  const unsigned int d = 3;
  Triangulation<d>   tr;

  // Create a pixel mask for 3d programmatically, rather than via
  // explicitly spelling it out as in the 2d case above:
  bool pixel_mask[7][5][9];
  for (unsigned int x = 0; x < 7; ++x)
    for (unsigned int y = 0; y < 5; ++y)
      for (unsigned int z = 0; z < 9; ++z)
        pixel_mask[x][y][z] =
          ((x % 2) == 1 && (y % 2 == 1) && (z % 2 == 1) ? false : true);

  GridGenerator::cheese(tr, Table<3, bool>(7, 5, 9, &pixel_mask[0][0][0]), 1.);

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
  dim_2a(logfile);
  dim_3(logfile);
}
