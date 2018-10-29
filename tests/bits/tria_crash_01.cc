// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2018 by the deal.II authors
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
