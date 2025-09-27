// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Check that we can create iterators filtered using ActiveFEIndexEqualTo.
//
// Compared to filtered_iterator_07, this test checks the availability
// of operator| to create the filter.


#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_generator.h>

#include <deal.II/hp/fe_collection.h>

#include "../tests.h"


template <int dim, int spacedim>
void
test()
{
  Triangulation<dim, spacedim> tria;
  GridGenerator::hyper_cube(tria, 0, 1, true);
  tria.refine_global(1);

  const FE_Q<dim, spacedim> fe(1);

  DoFHandler<dim, spacedim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  int n_cells_visited = 0;
  for (const auto &cell : dof_handler.active_cell_iterators() |
                            IteratorFilters::ActiveFEIndexEqualTo(0))
    {
      ++n_cells_visited;
      (void)cell;
    }

  Assert(n_cells_visited == Utilities::fixed_power<dim>(2),
         ExcMessage("Wrong number of cells visited."));

  n_cells_visited = 0;
  for (const auto &cell : dof_handler.active_cell_iterators() |
                            IteratorFilters::ActiveFEIndexEqualTo(1))
    {
      ++n_cells_visited;
      (void)cell;
    }

  Assert(n_cells_visited == 0, ExcMessage("Wrong number of cells visited."));
}


template <int dim, int spacedim>
void
test_hp()
{
  Triangulation<dim, spacedim> tria;
  GridGenerator::hyper_cube(tria, 0, 1, true);
  tria.refine_global(1);

  const hp::FECollection<dim, spacedim> fe_collection{FE_Q<dim, spacedim>(1),
                                                      FE_Q<dim, spacedim>(1)};

  DoFHandler<dim, spacedim> dof_handler(tria);

  for (auto &cell : dof_handler.active_cell_iterators())
    if (cell == dof_handler.begin_active())
      cell->set_active_fe_index(1);
    else
      cell->set_active_fe_index(0);

  dof_handler.distribute_dofs(fe_collection);

  int n_cells_visited = 0;
  for (const auto &cell : dof_handler.active_cell_iterators() |
                            IteratorFilters::ActiveFEIndexEqualTo(0))
    {
      ++n_cells_visited;
      (void)cell;
    }

  Assert(n_cells_visited == Utilities::fixed_power<dim>(2) - 1,
         ExcMessage("Wrong number of cells visited."));

  n_cells_visited = 0;
  for (const auto &cell : dof_handler.active_cell_iterators() |
                            IteratorFilters::ActiveFEIndexEqualTo(1))
    {
      ++n_cells_visited;
      (void)cell;
    }

  Assert(n_cells_visited == 1, ExcMessage("Wrong number of cells visited."));
}


int
main()
{
  initlog();

  test<2, 2>();
  test<2, 3>();
  test<3, 3>();

  test_hp<2, 2>();
  test_hp<2, 3>();
  test_hp<3, 3>();

  deallog << "OK" << std::endl;
}
