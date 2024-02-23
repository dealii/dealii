// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


/*
 * Compare find_active_cell_around_point between the default variant and the
 * one with a cache.
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
        // avoid roundoff issues of point at exactly the vertex
        p2[d] = (a + b) * (0.5 - 1e-8);
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
