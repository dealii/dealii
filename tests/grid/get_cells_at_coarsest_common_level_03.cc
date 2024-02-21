// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Test GridTools::get_cells_at_coarsest_common_level()



#include <deal.II/base/tensor.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"



template <int dim>
void
test()
{
  Triangulation<dim> triangulation(
    Triangulation<dim>::limit_level_difference_at_vertices);

  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global(1);
  // Adaptive refinement ... refine the first cell
  triangulation.begin_active()->set_refine_flag();
  triangulation.execute_coarsening_and_refinement();

  // now extract patches and print every fifth of them
  unsigned int index = 0;
  for (typename Triangulation<dim>::active_cell_iterator cell =
         triangulation.begin_active();
       cell != triangulation.end();
       ++cell, ++index)
    {
      std::vector<typename Triangulation<dim>::active_cell_iterator>
        patch_cells =
          GridTools::get_patch_around_cell<Triangulation<dim>>(cell);

      std::vector<typename Triangulation<dim>::cell_iterator> coarse_cells =
        GridTools::get_cells_at_coarsest_common_level<Triangulation<dim>>(
          patch_cells);

      deallog << "coarse_ cells " << cell << ": ";
      for (unsigned int i = 0; i < coarse_cells.size(); ++i)
        deallog << coarse_cells[i] << ' ';
      deallog << std::endl;
    }
}


int
main()
{
  initlog();

  deallog.push("1d");
  test<1>();
  deallog.pop();

  deallog.push("2d");
  test<2>();
  deallog.pop();

  deallog.push("3d");
  test<3>();
  deallog.pop();
}
