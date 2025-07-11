// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Create a Surface_mesh from each cell of an hyper ball, to check for
// orientation issues that may arise in Hexaedrons with non-standard face
// orientations.

#include <deal.II/base/config.h>

#include <deal.II/cgal/surface_mesh.h>
#include <deal.II/cgal/utilities.h>

#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"

using K         = CGAL::Exact_predicates_inexact_constructions_kernel;
using CGALPoint = CGAL::Point_3<K>;

int
main()
{
  initlog();
  Triangulation<3> tria;
  GridGenerator::hyper_ball(tria);

  for (const auto &cell : tria.active_cell_iterators())
    {
      CGAL::Surface_mesh<CGALPoint> mesh;
      CGALWrappers::dealii_cell_to_cgal_surface_mesh(
        cell, StaticMappingQ1<3>::mapping, mesh);
      deallog << "cell: " << cell << " is valid: " << (int)mesh.is_valid()
              << std::endl;
    }
}
