// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2017 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

// GridTools::regularize_corner_cells

#include "../tests.h"
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>

template <int dim>
void
test()
{
  const SphericalManifold<dim> m0;
  Triangulation<dim>           tria;
  GridGenerator::hyper_cube(tria, -1, 1);
  tria.set_all_manifold_ids_on_boundary(0);
  tria.set_manifold(0, m0);

  GridTools::regularize_corner_cells(tria);

  GridOut grid_out;
  grid_out.write_msh(tria, deallog.get_file_stream());
}

int
main()
{
  initlog();

  test<2>();

  return 0;
}
