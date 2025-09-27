// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2006 - 2021 by the deal.II authors
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

template <int dim>
void
test()
{
  Triangulation<dim> tria_1, tria_2, tria_3;
  GridGenerator::hyper_cube(tria_1);
  GridGenerator::hyper_cube(tria_2);

  // fill tria_3 with something, to
  // make sure that the function we
  // call later can deal with prior
  // content
  GridGenerator::hyper_cube(tria_3);

  // refine once, then refine first
  // cell
  tria_1.refine_global(1);
  tria_1.begin_active()->set_refine_flag();
  tria_1.execute_coarsening_and_refinement();

  // similar for second grid, but
  // different cell
  tria_2.refine_global(1);
  (std::next(tria_2.begin_active()))->set_refine_flag();
  tria_2.execute_coarsening_and_refinement();

  GridGenerator::create_union_triangulation(tria_1, tria_2, tria_3);

  GridOut().write_gnuplot(tria_3, deallog.get_file_stream());

  deallog << "     Total number of cells        = " << tria_3.n_cells()
          << std::endl
          << "     Total number of active cells = " << tria_3.n_active_cells()
          << std::endl;
}


int
main()
{
  initlog();
  deallog << std::setprecision(2);

  test<1>();
  test<2>();
  test<3>();

  return 0;
}
