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



#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"


// testcase=0:
// create two cubes; translate them so that no vertices overlap
//
// testcase=1:
// create two cubes; translate them so that a whole face overlaps
//
// testcase=2:
// create two cubes; translate them so that exactly one vertices overlaps
template <int dim>
void
test(const int testcase)
{
  Triangulation<dim> tria_1, tria_2, tria_3;
  GridGenerator::hyper_cube(tria_1);
  GridGenerator::hyper_cube(tria_2);
  Point<dim> shift;
  switch (testcase)
    {
      case 0:
        shift[0] = 2;
        break;
      case 1:
        shift[0] = 1;
        break;
      case 2:
        for (unsigned int d = 0; d < dim; ++d)
          shift[d] = 1;
        break;
      default:
        DEAL_II_NOT_IMPLEMENTED();
    }
  GridTools::shift(shift, tria_2);

  // fill tria_3 with something, to
  // make sure that the function we
  // call later can deal with prior
  // content
  GridGenerator::hyper_cube(tria_3);

  // now merge triangulations
  GridGenerator::merge_triangulations(tria_1, tria_2, tria_3);

  GridOut().write_gnuplot(tria_3, deallog.get_file_stream());

  deallog << "     Total number of cells        = " << tria_3.n_cells()
          << std::endl
          << "     Total number of vertices = " << tria_3.n_used_vertices()
          << std::endl;
}


int
main()
{
  initlog();
  deallog << std::setprecision(2);

  for (unsigned int t = 0; t < 3; ++t)
    {
      test<2>(t);
      test<3>(t);
    }

  return 0;
}
