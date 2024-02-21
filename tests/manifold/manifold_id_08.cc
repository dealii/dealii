// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Make sure that Triangulation::get_manifold_ids() work correctly

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>
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

  tria.set_manifold(1, FlatManifold<dim>());
  tria.set_manifold(2, FlatManifold<dim>());
  tria.set_manifold(3, FlatManifold<dim>());

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
