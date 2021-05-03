// ---------------------------------------------------------------------
//
// Copyright (C)  2015 - 2016  by the deal.II authors
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



// Test GridTools::build_triangulation_from_patch ()



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
