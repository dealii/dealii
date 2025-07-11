// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// check filtered iterators using multiple predicates

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <set>

#include "../tests.h"

using active_cell_iterator = Triangulation<2>::active_cell_iterator;

void
test()
{
  Triangulation<2> tria;
  GridGenerator::hyper_cube(tria, -1, 1);
  tria.refine_global(1);
  tria.begin_active()->set_refine_flag();
  tria.execute_coarsening_and_refinement();
  tria.refine_global(2);

  // we now have a number of cells,
  // flag them with some subdomain
  // ids based on their position, in
  // particular we take the quadrant
  // (octant)
  active_cell_iterator cell = tria.begin_active(), endc = tria.end();
  for (unsigned int i = 0; cell != endc; ++cell)
    {
      unsigned int subdomain = i % 3;

      cell->set_subdomain_id(subdomain);
      ++i;
    };

  // Count the cells that are on the boundary and have a subdomain_id of 0
  std::set<active_cell_iterator> cell_set;
  for (cell = tria.begin_active(); cell != endc; ++cell)
    if ((cell->subdomain_id() == 0) && (cell->at_boundary()))
      cell_set.insert(cell);


  unsigned int n_filtered_cells = 0;
  for (auto filtered_cell :
       filter_iterators(tria.active_cell_iterators(),
                        IteratorFilters::AtBoundary(),
                        IteratorFilters::SubdomainEqualTo(0)))
    {
      AssertThrow(cell_set.count(filtered_cell) == 1,
                  ExcMessage("Wrong cell filtered."));
      ++n_filtered_cells;
    }
  AssertThrow(n_filtered_cells == cell_set.size(),
              ExcMessage("Filtered cells missing."));
}

int
main()
{
  initlog();
  deallog << std::setprecision(4);

  test();

  deallog << "OK" << std::endl;

  return 0;
}
