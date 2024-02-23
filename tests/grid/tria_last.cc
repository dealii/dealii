// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2002 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// when we coarsen away a whole level, the level data structures
// remain but contains only unused cells. make sure that tria.last()
// and tria.last_active() still produce something sensible in that
// case

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"

template <int dim>
void
test()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(3);
  for (unsigned int i = 0; i < 2; ++i)
    {
      for (typename Triangulation<dim>::active_cell_iterator cell =
             tria.begin_active(2);
           cell != tria.end();
           ++cell)
        cell->set_coarsen_flag();
      tria.execute_coarsening_and_refinement();
    }

  deallog << tria.n_levels() << ' ' << tria.n_global_levels() << ' '
          << tria.last() << ' ' << tria.last_active() << std::endl;
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
