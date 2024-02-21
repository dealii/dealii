// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2001 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// GridTools::regularize_corner_cells on more complicated mesh

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>

#include "../tests.h"

int
main()
{
  initlog();

  Point<2> p0(0, 0), p1(4, 1);
  Point<2> c0(.1, .5), c1(3.9, .5);

  SphericalManifold<2> m0(c0);
  SphericalManifold<2> m1(c1);

  Triangulation<2>          tria;
  std::vector<unsigned int> subdivisions(2);
  subdivisions[0] = 4;
  subdivisions[1] = 1;
  GridGenerator::subdivided_hyper_rectangle(tria, subdivisions, p0, p1, true);

  GridTools::copy_boundary_to_manifold_id(tria);

  for (const auto bid : tria.get_boundary_ids())
    tria.set_manifold(bid, FlatManifold<2>());

  tria.set_manifold(0, m0);
  tria.set_manifold(1, m1);

  GridTools::regularize_corner_cells(tria);
  tria.refine_global(1);

  GridOut grid_out;
  grid_out.write_msh(tria, deallog.get_file_stream());

  return 0;
}
