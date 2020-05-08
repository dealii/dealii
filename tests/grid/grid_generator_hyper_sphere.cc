// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2019 by the deal.II authors
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
