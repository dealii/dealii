// ---------------------------------------------------------------------
//
// Copyright (C) 2022 by the deal.II authors
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

// Create a Surface_mesh from each cell of an hyper ball, to check for
// orientation issues that may arise in Hexaedrons with non-standard face
// orientations.

#include <deal.II/base/config.h>

#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/cgal/surface_mesh.h>
#include <deal.II/cgal/utilities.h>

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
