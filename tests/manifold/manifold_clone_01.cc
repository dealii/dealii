// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Check that now it is possible to build a mesh and a manifold in the
// "wrong" order.

#include <deal.II/base/utilities.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>

#include "../tests.h"

int
main()
{
  initlog();

  // We intentionally create the tria before the manifold, to make sure
  // that now it is possible to do so, since manifold is cloned internally
  // in the Triangulation

  Triangulation<2> tria;

  Point<2>                      center(1.0, 2.0);
  const SphericalManifold<2, 2> manifold(center);

  GridGenerator::hyper_ball(tria, center);
  GridTools::copy_boundary_to_manifold_id(tria);

  tria.set_manifold(0, manifold);
  tria.refine_global(1);

  GridOut go;
  go.write_msh(tria, deallog.get_file_stream());
}
