// ---------------------------------------------------------------------
//
// Copyright (C) 2019 by the deal.II authors
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

// Validate GridTools::assign_co_dimensional_manifold_indicators

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"


template <int dim, int spacedim>
void
test()
{
  deallog << "dim = " << dim << ", spacedim = " << spacedim << std::endl;

  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(1);
  tria.begin_active()->set_refine_flag();
  tria.execute_coarsening_and_refinement();

  tria.begin_active()->set_manifold_id(1);

  auto validate = [](const std::set<types::manifold_id> &set) {
    if (set.size() > 1)
      return 1;
    else
      return 0;
  };

  GridTools::assign_co_dimensional_manifold_indicators(tria, validate);
  for (auto &cell : tria.active_cell_iterators())
    {
      deallog << "Cell: " << (int)cell->manifold_id() << std::endl;
      if (dim > 1)
        for (unsigned int l = 0; l < GeometryInfo<dim>::lines_per_cell; ++l)
          deallog << "Line " << l << ", " << (int)cell->line(l)->manifold_id()
                  << std::endl;
      if (dim > 2)
        for (unsigned int l = 0; l < GeometryInfo<dim>::quads_per_cell; ++l)
          deallog << "Quad " << l << ", " << (int)cell->quad(l)->manifold_id()
                  << std::endl;
    }
}


int
main()
{
  initlog();

  test<1, 1>();
  test<1, 2>();
  test<2, 2>();
  test<2, 3>();
  test<3, 3>();

  return 0;
}
