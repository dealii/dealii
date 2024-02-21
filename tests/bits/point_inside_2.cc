// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2003 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check TriaAccessor<3>::point_inside for a cell that is _not_
// aligned with the coordinate axes
//
// this program is a modified version of one by Joerg Weimar,
// TU Braunschweig

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"


template <int dim>
void
check()
{
  // use a rather complicated way to
  // generate a new mesh with one
  // moved vertex -- one could move
  // it by hand in the original
  // triangulation, but this is how
  // the original testcase came in,
  // and it tests some more stuff, so
  // why change things
  Triangulation<dim> triangulation;

  // we generate a hyper_cube and
  // modify one vertex.
  GridGenerator::hyper_cube(triangulation);

  // Now get the cell
  const typename Triangulation<dim>::cell_iterator cell = triangulation.begin();
  cell->vertex(0)[0]                                    = -1.;

  // and test it.
  double testcoord[14][3] = {{0.5, 0.5, 0.5},
                             {2, 0.5, 0.5},
                             {0.5, 2, 0.5},
                             {0.5, 0.5, 2},
                             {-2, 0.5, 0.5},
                             {0.5, -2, 0.5},
                             {0.5, 0.5, -2},
                             {0.9, 0.9, 0.9},
                             {1.0, 0.5, 0.5},
                             {0.9999999, 0.5, 0.5},
                             {1.0000001, 0.5, 0.5},
                             {-0.1, 0.1, 0.1},
                             {-0.24, 0.5, 0.5},
                             {-0.26, 0.5, 0.5}};

  const bool  expected2d[] = {1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1};
  const bool  expected3d[] = {1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0};
  const bool *expected     = dim == 2 ? expected2d : expected3d;
  for (int i = 0; i < 14; ++i)
    {
      Point<dim> testpoint;
      testpoint[0] = testcoord[i][0];
      testpoint[1] = testcoord[i][1];
      if (dim == 3)
        testpoint[2] = testcoord[i][2];

      bool res = cell->point_inside(testpoint);
      deallog << testpoint << "  \t inside " << res << " expected "
              << expected[i] << std::endl;
      Assert(res == expected[i], ExcInternalError());
    }
}


int
main()
{
  initlog();

  check<2>();
  check<3>();
}
