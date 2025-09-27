// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2014 - 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Check range-based for loops for triangulations

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <string>

#include "../tests.h"


template <int dim>
void
check()
{
  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr);
  tr.refine_global(2);


  {
    // set flags on active cells
    tr.clear_user_flags();
    for (auto &cell : tr.active_cell_iterators())
      cell->set_user_flag();

    // now verify that it is really only the active cells
    for (auto &cell : tr.cell_iterators())
      AssertThrow(cell->user_flag_set() == !cell->has_children(),
                  ExcInternalError());
  }

  // now do the same again for all levels of the triangulation
  for (unsigned int l = 0; l < tr.n_levels(); ++l)
    {
      tr.clear_user_flags();
      for (auto &cell : tr.active_cell_iterators_on_level(l))
        cell->set_user_flag();

      for (auto &cell : tr.cell_iterators_on_level(l))
        AssertThrow(cell->user_flag_set() == !cell->has_children(),
                    ExcInternalError());

      for (auto &cell : tr.cell_iterators())
        AssertThrow((cell->user_flag_set() == !cell->has_children()) ||
                      (l != (unsigned int)cell->level()),
                    ExcInternalError());
    }

  deallog << "OK" << std::endl;
}


int
main()
{
  deal_II_exceptions::disable_abort_on_exception();

  initlog();

  check<2>();
  check<3>();
}
