// ---------------------------------------------------------------------
//
// Copyright (C) 2014 by the deal.II authors
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



// Check range-based for loops for triangulations

#include "../tests.h"
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/base/logstream.h>

#include <fstream>
#include <string>


template<int dim>
void check()
{
  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr);
  tr.refine_global(2);


  {
    // set flags on active cells
    tr.clear_user_flags ();
    for (auto cell : tr.active_cell_iterators())
      cell->set_user_flag();

    // now verify that it is really only the active cells
    for (auto cell : tr.cell_iterators())
      Assert (cell->user_flag_set() == !cell->has_children(),
	      ExcInternalError());
  }

  // now do the same again for all levels of the triangulation
  for (unsigned int l=0; l<tr.n_levels(); ++l)
    {
      tr.clear_user_flags ();
      for (auto cell : tr.active_cell_iterators_on_level(l))
	cell->set_user_flag();

      for (auto cell : tr.cell_iterators_on_level(l))
	Assert (cell->user_flag_set() == !cell->has_children(),
		ExcInternalError());
      
      for (auto cell : tr.cell_iterators())
	Assert ((cell->user_flag_set() == !cell->has_children())
		||
		(l != cell->level()),
		ExcInternalError());
    }

  deallog << "OK" << std::endl;
}


int main()
{
  deal_II_exceptions::disable_abort_on_exception();

  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  check<2>();
  check<3>();
}
