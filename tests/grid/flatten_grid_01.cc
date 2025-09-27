// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2014 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Generate a grid, refine it once, flatten it and output the result.

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"

template <int dim, int spacedim1, int spacedim2>
void
test()
{
  deallog << "Testing <" << dim << ',' << spacedim1 << "> VS <" << dim << ','
          << spacedim2 << '>' << std::endl;

  Triangulation<dim, spacedim1> tria1;
  GridGenerator::hyper_cube(tria1);
  tria1.refine_global(1);

  Triangulation<dim, spacedim2> tria2;
  GridGenerator::flatten_triangulation(tria1, tria2);
  GridOut go;
  go.write_msh(tria2, deallog.get_file_stream());
}

int
main()
{
  initlog();
  test<1, 1, 1>();
  test<1, 1, 2>();
  test<1, 1, 3>();
  //
  test<1, 2, 1>();
  test<1, 2, 2>();
  test<1, 2, 3>();
  //
  test<1, 3, 1>();
  test<1, 3, 2>();
  test<1, 3, 3>();
  //
  test<2, 2, 2>();
  test<2, 2, 3>();
  //
  test<2, 3, 2>();
  test<2, 3, 3>();
  //
  test<3, 3, 3>();
}
