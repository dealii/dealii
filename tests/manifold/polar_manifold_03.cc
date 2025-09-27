// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test the push_forward and pull_back mechanisms

#include "../tests.h"


// all include files you need here
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>


// Helper function
template <int dim, int spacedim>
void
test(unsigned int ref = 1)
{
  deallog << "Testing dim " << dim << ", spacedim " << spacedim << std::endl;

  PolarManifold<dim, spacedim> manifold;

  Triangulation<dim, spacedim> tria;
  Point<spacedim>              p0;
  Point<spacedim>              p1;
  p0[0] = .2;
  p1[0] = 1;
  p0[1] = .1;

  if (spacedim == 2)
    {
      p1[1] = 2 * numbers::PI - .1; // theta
    }
  else if (spacedim == 3)
    {
      p1[1] = numbers::PI - .1;
      p1[2] = 2 * numbers::PI - .1;
    }

  GridGenerator::hyper_rectangle(tria, p0, p1);
  tria.refine_global(3);

  const std::vector<Point<spacedim>> &vertices = tria.get_vertices();

  for (unsigned int i = 0; i < vertices.size(); ++i)
    {
      Point<spacedim> p0 = manifold.push_forward(vertices[i]);
      Point<spacedim> p1 = manifold.pull_back(p0);

      if (p1.distance(vertices[i]) > 1e-10)
        deallog << "ERROR! d: " << p1.distance(vertices[i]) << " - " << p1
                << " != " << vertices[i] << std::endl;
    }
}

int
main()
{
  initlog();

  test<2, 2>();
  test<3, 3>();

  return 0;
}
