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
  triangulation.refine_global(2);


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

      // in uniform refinemenet(without any hanging nodes), we expect the size of vector of patch_cells for each cell be
      // equal to the vector of coarse_cells around that cell
      Assert(patch_cells.size() == coarse_cells.size(), ExcInternalError());

      // we sort both vectors to ensure that we are comparing the correct cells to each other
      std::sort(patch_cells.begin(), patch_cells.end());
      std::sort(coarse_cells.begin(), coarse_cells.end());

      for (unsigned int i = 0; i < patch_cells.size(); ++i)
        Assert(patch_cells[i] == coarse_cells[i], ExcInternalError());

      deallog << "coarse_ cells " << cell << ": ";
      for (unsigned int i = 0; i < coarse_cells.size(); ++i)
        deallog << coarse_cells[i] << ' ';
      deallog << std::endl;
    }

  deallog << "OKay" << std::endl;
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
