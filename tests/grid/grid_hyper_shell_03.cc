// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2007 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// generate and refine a hyper shell in 3d with 12 cells. The other tests of
// this series test the other possibilities for the number of cells

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
  static const SphericalManifold<dim> boundary(center);
  tria.set_manifold(0, boundary);
  tria.refine_global(1);

  GridOut grid_out;
  grid_out.write_gnuplot(tria, deallog.get_file_stream());
}


int
main()
{
  initlog();
  deallog << std::setprecision(3);

  check<3>(.5, 1, 12);
}
