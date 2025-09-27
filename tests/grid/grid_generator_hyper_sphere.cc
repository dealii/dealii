// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Test GridGenerator::hyper_sphere()

#include <deal.II/base/tensor.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"



template <int dim, int spacedim>
void
test(std::ostream &out)
{
  Triangulation<dim, spacedim> triangulation;

  static SphericalManifold<dim, spacedim> surface_description;

  GridGenerator::hyper_sphere(triangulation);

  triangulation.set_all_manifold_ids(0);
  triangulation.set_manifold(0, surface_description);
  triangulation.refine_global(3);

  GridOut go;
  go.write_gnuplot(triangulation, out);
}


int
main()
{
  initlog();

  test<1, 2>(deallog.get_file_stream());
  test<2, 3>(deallog.get_file_stream());
}
