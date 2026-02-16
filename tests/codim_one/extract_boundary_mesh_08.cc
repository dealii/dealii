// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2010 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


// A test that used to fail because of mis-oriented faces

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"



void
test()
{
  const int dim = 3;

  Triangulation<dim> triangulation;
  GridGenerator::cylinder(triangulation);

  GridOut().write_gnuplot(triangulation, deallog.get_file_stream());

  Triangulation<dim - 1, dim> triangulation_surface;

  for (const auto bid : triangulation.get_manifold_ids())
    if (bid != numbers::flat_manifold_id)
      triangulation_surface.set_manifold(bid, FlatManifold<2, 3>());

  GridGenerator::extract_boundary_mesh(triangulation, triangulation_surface);
  triangulation_surface.refine_global(2);

  GridOut().write_gnuplot(triangulation_surface, deallog.get_file_stream());

  deallog << triangulation_surface.n_used_vertices() << std::endl;
  deallog << triangulation_surface.n_active_cells() << std::endl;
}


int
main()
{
  initlog();

  test();
}
