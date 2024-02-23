// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2007 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// GridGenerator::hyper_shell colorized the faces but forgot the
// edges. This is not useful because the colorization is usually done
// so that one can attach a boundary or manifold object to these parts
// of the boundary

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
  GridGenerator::hyper_shell(tria, center, r1, r2, n, true);
  static const SphericalManifold<dim> boundary(center);
  tria.set_manifold(0, boundary);

  for (typename Triangulation<dim>::cell_iterator cell = tria.begin();
       cell != tria.end();
       ++cell)
    for (const unsigned int f : GeometryInfo<dim>::face_indices())
      if (cell->face(f)->at_boundary())
        for (unsigned int l = 0; l < GeometryInfo<dim>::lines_per_face; ++l)
          AssertThrow(cell->face(f)->line(l)->boundary_id() ==
                        cell->face(f)->boundary_id(),
                      ExcInternalError());

  deallog << "OK" << std::endl;
}


int
main()
{
  initlog();
  deallog << std::setprecision(3);

  check<3>(.5, 1, 6);
  check<3>(.5, 1, 12);
  check<3>(.5, 1, 96);
}
