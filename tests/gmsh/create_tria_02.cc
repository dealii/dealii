// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Create a hyper ball, refine it, extract an iges of the boundary,
// and create a new mesh using gmsh from that iges file.

#include <deal.II/gmsh/utilities.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include <deal.II/opencascade/manifold_lib.h>
#include <deal.II/opencascade/utilities.h>

#include "../tests.h"

int
main()
{
  initlog();

  Triangulation<2> tria(Triangulation<2>::MeshSmoothing::none,
                        /*check_for_distorted_cells*/ true);
  GridGenerator::hyper_ball(tria);

  tria.refine_global(4);

  auto curves = OpenCASCADE::create_curves_from_triangulation_boundary(tria);
  tria.clear();

  Gmsh::create_triangulation_from_boundary_curve(curves[0], tria);

  // The grid created depends on the OpenCASCADE and Gmsh version used.
  // Hence, only check that the resulting mesh is non-empty and the cells
  // are not distorted.
  AssertThrow(tria.n_cells() > 0, ExcInternalError());

  deallog << "OK" << std::endl;

  return 0;
}
