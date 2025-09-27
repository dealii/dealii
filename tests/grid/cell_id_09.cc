// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


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
