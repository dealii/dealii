// ---------------------------------------------------------------------
//
// Copyright (C) 2022 by the deal.II authors
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



// check CellAccessor::is_locally_owned_on_level and is_ghost_on_level for a
// serial triangulation (where trivially all cells are locally owned and none
// is ghosted/artificial), going through a code path that is not tested
// otherwise

#include <deal.II/base/utilities.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>

#include "../tests.h"



template <int dim>
void
test()
{
  Triangulation<dim> tr;

  GridGenerator::subdivided_hyper_cube(tr, 2);

  for (const auto &cell : tr.cell_iterators())
    {
      if (cell->is_locally_owned_on_level())
        {
          deallog << cell << ": locally owned" << std::endl;
          Assert(!cell->is_ghost_on_level() && !cell->is_artificial_on_level(),
                 ExcInternalError());
        }
      if (cell->is_ghost_on_level())
        {
          deallog << cell << ": ghost" << std::endl;
          Assert(!cell->is_locally_owned_on_level() &&
                   !cell->is_artificial_on_level(),
                 ExcInternalError());
        }
      if (cell->is_artificial_on_level())
        {
          deallog << cell << ": artificial" << std::endl;
          Assert(!cell->is_locally_owned_on_level() &&
                   !cell->is_ghost_on_level(),
                 ExcInternalError());
        }
    }
}


int
main()
{
  initlog();
  test<2>();
}
