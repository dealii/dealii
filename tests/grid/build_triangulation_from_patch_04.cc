// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Test GridTools::build_triangulation_from_patch () with a distorted
// triangulation to make sure the vertices are captured properly when the patch
// triangulation does not come from the standard deal.II refinement strategy.
// This mesh is not uniform and so the local_triangulations return will have
// additional cells that are not in the desired patch but that make sure we have
// a valid triangulation


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
  triangulation.begin_active()->set_refine_flag();
  triangulation.execute_coarsening_and_refinement();

  // Apply distortion here, note that distort_random
  // is designed to be reproducible, each time you apply
  // it to the same triangulation, you get the same
  // random distortion.  Thus the output is consistent
  // for the test to be run multiple times.
  GridTools::distort_random(0.1, triangulation, false);


  unsigned int index = 0;
  for (typename Triangulation<dim>::active_cell_iterator cell =
         triangulation.begin_active();
       cell != triangulation.end();
       ++cell, ++index)
    {
      std::vector<typename Triangulation<dim>::active_cell_iterator>
        patch_cells =
          GridTools::get_patch_around_cell<Triangulation<dim>>(cell);

      Triangulation<dim> local_triangulation(
        Triangulation<dim>::limit_level_difference_at_vertices);
      std::map<typename Triangulation<dim>::active_cell_iterator,
               typename Triangulation<dim>::active_cell_iterator>
        patch_to_global_tria_map;

      GridTools::build_triangulation_from_patch<Triangulation<dim>>(
        patch_cells, local_triangulation, patch_to_global_tria_map);

      deallog << "patch_cells " << cell << ": ";
      for (unsigned int i = 0; i < patch_cells.size(); ++i)
        deallog << patch_cells[i] << ' ';
      deallog << std::endl;

      deallog << "local_triangulation " << cell << ": ";
      for (const auto &tria_cell : local_triangulation.active_cell_iterators())
        {
          deallog << "   " << tria_cell << " user flag check:  "
                  << (tria_cell->user_flag_set() ? " (+) " : " (-) ")
                  << std::endl;
          for (const unsigned int v : GeometryInfo<dim>::vertex_indices())
            {
              deallog << "  vertices for cell  " << tria_cell << " : "
                      << tria_cell->vertex(v) << std::endl;
            }
        }
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
