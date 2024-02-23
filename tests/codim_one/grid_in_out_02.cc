// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2008 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// like grid_in_out but write in gnuplot format

#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include <string>

#include "../tests.h"

template <int dim, int spacedim>
void
test(std::string filename)
{
  Triangulation<dim, spacedim> tria;
  GridIn<dim, spacedim>        gi;
  gi.attach_triangulation(tria);
  std::ifstream in(filename);
  gi.read_ucd(in);

  GridOut grid_out;
  grid_out.write_gnuplot(tria, deallog.get_file_stream());
}

int
main()
{
  initlog();

  test<2, 3>(SOURCE_DIR "/grids/square.inp");
  test<2, 3>(SOURCE_DIR "/grids/sphere_1.inp");

  return 0;
}
