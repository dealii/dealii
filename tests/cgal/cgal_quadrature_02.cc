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

// Intersect reference cells and compute the volume of the intersection by
// integration.

#include <deal.II/base/config.h>

#include <deal.II/cgal/utilities.h>

#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <CGAL/IO/File_medit.h>
#include <CGAL/IO/io.h>
#include <CGAL/Polygon_mesh_processing/measure.h>

#include "../tests.h"

// Compute volumes of intersections of polyhedra by splitting each intersection
// into simplices.

using namespace CGALWrappers;

void
test()
{
  using namespace ReferenceCells;
  std::vector<std::pair<ReferenceCell, ReferenceCell>> ref_pairs = {
    {Hexahedron, Pyramid},
    {Hexahedron, Tetrahedron},
    {Tetrahedron, Pyramid},
    {Pyramid, Pyramid}};

  constexpr int  degree  = 3;
  constexpr auto bool_op = BooleanOperation::compute_intersection;
  for (const auto &pair : ref_pairs)
    {
      const auto       ref_cell0 = pair.first;
      const auto       ref_cell1 = pair.second;
      Triangulation<3> tria0;
      Triangulation<3> tria1;
      GridGenerator::reference_cell(tria0, ref_cell0);
      GridGenerator::reference_cell(tria1, ref_cell1);
      const auto mapping0 = ref_cell0.template get_default_mapping<3>(1);
      const auto mapping1 = ref_cell1.template get_default_mapping<3>(1);
      const auto cell0    = tria0.begin_active();
      const auto cell1    = tria1.begin_active();

      // Shift each vertex by the same amount
      const unsigned int n_vertices = cell1->n_vertices();
      for (unsigned int i = 0; i < n_vertices; ++i)
        {
          cell1->vertex(i) += Point<3>(0.1, 0.1, 0.1);
        }

      auto test_quad = compute_quadrature_on_boolean_operation<3, 3, 3>(
        cell0, cell1, degree, bool_op, *mapping0, *mapping1);
      deallog << "Volume of poly with Quadrature: " << std::setprecision(12)
              << std::accumulate(test_quad.get_weights().begin(),
                                 test_quad.get_weights().end(),
                                 0.)
              << std::endl;
    }
}

int
main()
{
  initlog();
  test();
}
