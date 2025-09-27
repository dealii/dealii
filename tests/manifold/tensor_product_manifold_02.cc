// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test TensorProductManifold by refining and generating normals for
// a manually constructed cylinder hull.

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tensor_product_manifold.h>

#include "../tests.h"


void
test()
{
  std::ostream &out = deallog.get_file_stream();

  FunctionManifold<1, 1> F("x", "x");
  PolarManifold<2, 2>    G;

  TensorProductManifold<2, 1, 1, 1, 2, 2, 2> manifold(F, G);

  // make a hull of a cylinder
  Triangulation<2, 3> tria;
  {
    Triangulation<3, 3> volume_tria;
    GridGenerator::cylinder(volume_tria);
    const std::set<types::boundary_id> boundary_ids = {0};
    GridGenerator::extract_boundary_mesh(volume_tria, tria, boundary_ids);
  }
  tria.set_all_manifold_ids(0);
  tria.set_manifold(0, manifold);

  tria.refine_global(1);

  out << "set view equal xyz" << std::endl
      << "splot '-' with lines, '-' with vectors " << std::endl;
  GridOut().write_gnuplot(tria, out);
  out << 'e' << std::endl;

  Triangulation<2, 3>::active_cell_iterator it = tria.begin_active();
  for (; it != tria.end(); ++it)
    {
      Point<3>     p  = it->center(true);
      Tensor<1, 3> t1 = manifold.get_tangent_vector(p, it->vertex(0));
      Tensor<1, 3> t2 = manifold.get_tangent_vector(p, it->vertex(1));
      Tensor<1, 3> n  = cross_product_3d(t1, t2);
      n /= -n.norm();
      out << it->center() << ' ' << n << std::endl;
    }
  out << 'e' << std::endl;
}



int
main()
{
  initlog();

  test();

  return 0;
}
