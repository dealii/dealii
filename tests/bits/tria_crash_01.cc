// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2006 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// a test that checks for a crash introduced in the triangulation class in the
// last few days when fixing refine_and_coarsen_3d


#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>

#include "../tests.h"



bool
predicate(const Point<3> &p, const double diameter)
{
  return ((p[0] - .2) * (p[0] - .2) + (p[2] - p[1] / 4) * (p[2] - p[1] / 4) <
          diameter * diameter);
}


int
main()
{
  initlog();

  const unsigned int dim = 3;
  Triangulation<dim> tria;
  GridGenerator::cylinder(tria, 1, .7);
  tria.reset_all_manifolds();

  deallog << "n_cells=" << tria.n_active_cells() << std::endl;

  tria.refine_global(2);

  deallog << "n_cells=" << tria.n_active_cells() << std::endl;

  // build up a map of vertex indices
  // of boundary vertices to the new
  // boundary points
  std::map<unsigned int, Point<dim>> new_points;

  Triangulation<dim>::active_cell_iterator cell = tria.begin_active(),
                                           endc = tria.end();

  for (cell = tria.begin_active(); cell != endc; ++cell)
    if (predicate(cell->center(), cell->diameter()))
      cell->set_refine_flag();
  tria.execute_coarsening_and_refinement();

  deallog << "n_cells=" << tria.n_active_cells() << std::endl;


  for (cell = tria.begin_active(); cell != endc; ++cell)
    if (!predicate(cell->center(), cell->diameter()))
      cell->set_coarsen_flag();

  // make sure there really are no refinement
  // flags set
  tria.prepare_coarsening_and_refinement();
  for (cell = tria.begin_active(); cell != endc; ++cell)
    AssertThrow(!cell->refine_flag_set(), ExcInternalError());

  tria.execute_coarsening_and_refinement();

  deallog << "n_cells=" << tria.n_active_cells() << std::endl;
}
