// ---------------------------------------------------------------------
//
// Copyright (C) 2017 by the deal.II authors
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


/*
 * Given the number of refinements and the number of random points
 * it benchmarks the time needed to run the function FCT
 * which can be point_locator_D2 (or point_locator when it shall be written)
 */
#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_tools_cache.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"

template <int dim, int spacedim>
void
test()
{
  deallog << "Testing " << dim << ", " << spacedim << std::endl;
  // The hypercube is [a,b]^spacedim
  const double a = -0.3;
  const double b = 0.7;


  Triangulation<dim, spacedim> triangulation;
  GridGenerator::hyper_cube(triangulation, a, b);
  triangulation.refine_global(dim == 2 ? 2 : 1);

  MappingQ<dim> mapping(1);

  DoFHandler<spacedim> dof_handler(triangulation);

  GridTools::Cache<dim, spacedim> cache(triangulation, mapping);
  Point<spacedim>                 p1;
  Point<spacedim>                 p2;
  Point<spacedim>                 p3;
  for (unsigned int d = 0; d < spacedim; ++d)
    {
      p1[d] = a;
      if (d == 1)
        p2[d] = b;
      else
        p2[d] = (a + b) * 0.5;
      p3[d] = b;
    }

  std::pair<typename Triangulation<dim, spacedim>::active_cell_iterator,
            Point<dim>>
    res_1 = GridTools::find_active_cell_around_point(mapping, dof_handler, p1);
  std::pair<typename Triangulation<dim, spacedim>::active_cell_iterator,
            Point<dim>>
    res_2 = GridTools::find_active_cell_around_point(mapping, dof_handler, p2);
  std::pair<typename Triangulation<dim, spacedim>::active_cell_iterator,
            Point<dim>>
    res_3 = GridTools::find_active_cell_around_point(mapping, dof_handler, p3);

  deallog << "Standard find active cell around point: OK " << std::endl;

  // Failing with cache
  auto res_c_1 = GridTools::find_active_cell_around_point(cache, p1);
  auto res_c_2 = GridTools::find_active_cell_around_point(cache, p2);
  auto res_c_3 = GridTools::find_active_cell_around_point(cache, p3);

  deallog << "Cache find active cell around point: OK " << std::endl;

  if (res_1.first != res_c_1.first)
    deallog << "Different cells were found for p1" << std::endl;
  else if (res_2.first != res_c_2.first)
    deallog << "Different cells were found for p2" << std::endl;
  else if (res_3.first != res_c_3.first)
    deallog << "Different cells were found for p3" << std::endl;

  deallog << "Test: OK " << std::endl;
}

int
main()
{
  initlog();
  test<1, 1>();
  test<2, 2>();
  test<3, 3>();
}
