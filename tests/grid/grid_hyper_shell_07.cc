// ---------------------------------------------------------------------
//
// Copyright (C) 2020 by the deal.II authors
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


// test the various cell variants of GridGenerator::hyper_shell in 3d

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <iostream>

#include "../tests.h"

template <int dim>
void
check(double r1, double r2, unsigned int n)
{
  Point<dim>         center;
  Triangulation<dim> tria(Triangulation<dim>::none);
  GridGenerator::hyper_shell(tria, center, r1, r2, n);

  deallog << "Number of cells: " << tria.n_cells() << std::endl;

  unsigned int n_r1           = 0;
  unsigned int n_r2           = 0;
  unsigned int n_other_radius = 0;

  // ensure that all vertices of the mesh are either at r1 or at r2
  for (const auto v : tria.get_vertices())
    if (std::abs(v.norm_square() - r1 * r1) < 1e-12)
      ++n_r1;
    else if (std::abs(v.norm_square() - r2 * r2) < 1e-12)
      ++n_r2;
    else
      ++n_other_radius;

  deallog << "Number of vertices at inner radius: " << n_r1 << std::endl;
  deallog << "Number of vertices at outer radius: " << n_r2 << std::endl;
  deallog << "Number of vertices at other radius: " << n_other_radius
          << std::endl;
  deallog << std::endl;
}


int
main()
{
  initlog();

  check<3>(.8, 1, 6);
  check<3>(.8, 1, 12);
  check<3>(.8, 1, 24);
  check<3>(.8, 1, 48);
  check<3>(.8, 1, 96);
  check<3>(.8, 1, 192);
  check<3>(.8, 1, 384);
  check<3>(.8, 1, 768);
  check<3>(.8, 1, 1536);
  check<3>(.8, 1, 3072);
  check<3>(.8, 1, 6144);
}
