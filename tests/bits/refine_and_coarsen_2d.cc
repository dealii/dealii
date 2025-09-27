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



// check that if we take an locally refined mesh, refine it globally once,
// then coarsen it globally again, that we get the same mesh
//
// Triangulation::fix_coarsen_flags used to be too conservative in allowing
// cells to be coarsened (see today's changes.html entry)

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"



template <int dim>
void
check()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(2);
  tria.begin_active()->set_refine_flag();
  tria.execute_coarsening_and_refinement();

  // store which cells we have here
  std::vector<typename Triangulation<dim>::active_cell_iterator> cells;
  for (typename Triangulation<dim>::active_cell_iterator cell =
         tria.begin_active();
       cell != tria.end();
       ++cell)
    cells.push_back(cell);

  const unsigned int n_cells = tria.n_active_cells();
  deallog << n_cells << std::endl;

  // refine the mesh globally, then coarsen
  // it again globally
  tria.refine_global(1);

  for (typename Triangulation<dim>::active_cell_iterator cell =
         tria.begin_active();
       cell != tria.end();
       ++cell)
    cell->set_coarsen_flag();
  tria.execute_coarsening_and_refinement();


  // verify that we get the same cells again
  deallog << n_cells << ' ' << tria.n_active_cells() << std::endl;

  Assert(tria.n_active_cells() == n_cells, ExcInternalError());

  unsigned int index = 0;
  for (typename Triangulation<dim>::active_cell_iterator cell =
         tria.begin_active();
       cell != tria.end();
       ++cell, ++index)
    AssertThrow(cells[index] == cell, ExcInternalError());
}


int
main()
{
  initlog();

  check<2>();
}
