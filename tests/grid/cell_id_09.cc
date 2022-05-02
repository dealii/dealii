// ---------------------------------------------------------------------
//
// Copyright (C) 2020 by the deal.II authors
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


// verify getter functions of CellId

#include <deal.II/grid/cell_id.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"


template <int dim>
void
test(const unsigned int n_global_refinements)
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(n_global_refinements);

  // create a new CellId object from an existing one using just the getter
  // functions and verify both original and copy are equal
  for (const auto &cell : tria.active_cell_iterators())
    {
      CellId cell_id = cell->id();

      const auto coarse_cell_id = cell_id.get_coarse_cell_id();
      const auto child_indices  = cell_id.get_child_indices();

      CellId copy_id(coarse_cell_id,
                     child_indices.size(),
                     child_indices.data());

      AssertThrow(cell_id == copy_id, ExcInternalError());
    }

  deallog << "OK" << std::endl;
}


int
main()
{
  initlog();

  test<2>(3);
}
