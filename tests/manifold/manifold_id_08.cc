// ---------------------------------------------------------------------
//
// Copyright (C) 2021 by the deal.II authors
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

// Make sure that Triangulation::get_manifold_ids() work correctly

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"

template <int dim, int spacedim>
void
test(unsigned int ref = 1)
{
  deallog << "Testing dim=" << dim << ", spacedim=" << spacedim << std::endl;

  Triangulation<dim, spacedim> tria;
  GridGenerator::hyper_cube(tria);

  typename Triangulation<dim, spacedim>::active_cell_iterator cell;

  tria.begin_active()->set_manifold_id(1);
  for (const unsigned int f : tria.begin_active()->face_indices())
    tria.begin_active()->face(f)->set_manifold_id(2);

  tria.begin_active()->line(0)->set_manifold_id(3);

  tria.refine_global(ref);

  for (cell = tria.begin_active(); cell != tria.end(); ++cell)
    {
      deallog << "C: " << cell << ", manifold_id: " << (int)cell->manifold_id()
              << std::endl;
      for (const unsigned int f : cell->face_indices())
        deallog << "f: " << cell->face(f)
                << ", manifold_id: " << (int)cell->face(f)->manifold_id()
                << std::endl;
      if (dim == 3)
        for (const unsigned int l : cell->line_indices())
          deallog << "line: " << cell->line(l)
                  << ", manifold_id: " << (int)cell->line(l)->manifold_id()
                  << std::endl;
    }

  const auto manifold_ids = tria.get_manifold_ids();
  deallog << "All manifold ids: ";
  std::string sep = "";
  for (const auto &mid : manifold_ids)
    {
      deallog << sep << (int)mid;
      sep = ", ";
    }
  deallog << std::endl;
}

int
main()
{
  initlog();

  test<3, 3>(0);
  test<3, 3>(1);
}
