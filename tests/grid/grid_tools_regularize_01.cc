// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2001 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// GridTools::regularize_corner_cells

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>

#include "../tests.h"

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
