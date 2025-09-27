// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test that the flat manifold does what it should on a sphere surface.

#include "../tests.h"



// all include files you need here
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

// Helper function
template <int dim, int spacedim>
void
test(unsigned int ref = 1)
{
  PolarManifold<dim, spacedim> manifold;

  Triangulation<spacedim, spacedim> volume_tria;
  Triangulation<dim, spacedim>      tria;
  GridGenerator::hyper_ball(volume_tria);
  GridGenerator::extract_boundary_mesh(volume_tria, tria);

  typename Triangulation<dim, spacedim>::active_cell_iterator cell;

  for (cell = tria.begin_active(); cell != tria.end(); ++cell)
    cell->set_all_manifold_ids(1);

  for (cell = tria.begin_active(); cell != tria.end(); ++cell)
    {
      if ((fabs(cell->center()[0]) < 1e-10) &&
          (fabs(cell->center()[1]) < 1e-10))
        cell->set_all_manifold_ids(0);
    }

  tria.set_manifold(0, FlatManifold<dim, spacedim>());
  tria.set_manifold(1, manifold);
  tria.refine_global(2);

  GridOut gridout;
  gridout.write_msh(tria, deallog.get_file_stream());
}

int
main()
{
  initlog();

  test<1, 2>();
  test<2, 3>();

  return 0;
}
