// ---------------------------------------------------------------------
//
// Copyright (C) 2016 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


// check filtered iterators using multiple predicate

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/filtered_iterator.h>
#include <set>

typedef Triangulation<2>::active_cell_iterator active_cell_iterator;

void test ()
{
  Triangulation<2> tria;
  GridGenerator::hyper_cube(tria, -1, 1);
  tria.refine_global (1);
  tria.begin_active()->set_refine_flag ();
  tria.execute_coarsening_and_refinement ();
  tria.refine_global (2);

  // we now have a number of cells,
  // flag them with some subdomain
  // ids based on their position, in
  // particular we take the quadrant
  // (octant)
  active_cell_iterator cell = tria.begin_active (),
                       endc = tria.end ();
  for (unsigned int i=0; cell!=endc; ++cell)
    {
      unsigned int subdomain = i%3;

      cell->set_subdomain_id (subdomain);
      ++i;
    };

  // Count the cells that are on the boundary and have a subdomain_id of 0
  std::set<active_cell_iterator> cell_set;
  for (cell = tria.begin_active(); cell!=endc; ++cell)
    if ((cell->subdomain_id()==0) && (cell->at_boundary()))
      cell_set.insert(cell);


  unsigned int n_filtered_cells = 0;
  for (auto filtered_cell : filter_iterators(tria.active_cell_iterators(),
                                             IteratorFilters::AtBoundary(),
                                             IteratorFilters::SubdomainEqualTo(0)))
    {
      AssertThrow(cell_set.count(filtered_cell)==1, ExcMessage("Wrong cell filtered."));
      ++n_filtered_cells;
    }
  AssertThrow(n_filtered_cells==cell_set.size(), ExcMessage("Filtered cells missing."));
}

int main ()
{
  initlog();
  deallog << std::setprecision(4);
  deallog.threshold_double(1.e-10);

  test ();

  deallog << "OK" << std::endl;;

  return 0;
}

