/* ------------------------------------------------------------------------
 *
 * SPDX-License-Identifier: LGPL-2.1-or-later
 * Copyright (C) 2024 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * Part of the source code is dual licensed under Apache-2.0 WITH
 * LLVM-exception OR LGPL-2.1-or-later. Detailed license information
 * governing the source code and code contributions can be found in
 * LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
 *
 * ------------------------------------------------------------------------
 */

// Test that we get correctly oriented quadrilateral faces with mixed meshes.
// This verifies that we no longer assume that quadrilateral elements are always
// consistently oriented.

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_simplex_p.h>

#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/fe_values.h>
#include <deal.II/hp/mapping_collection.h>
#include <deal.II/hp/q_collection.h>

#include "../tests.h"

int
main()
{
  initlog();

  const std::vector<Point<2>> vertices{
    {0.0, 0.0}, {1.0, 0.0}, {0.0, 2.0}, {3.0, 1.0}, {3.0, 3.0}, {2.0, 3.0}};

  std::vector<CellData<2>> cells(3);
  cells[0].vertices = {0u, 1u, 2u};
  // 3, 5, 4 puts all faces in the standard orientation
  // 3, 4, 5 produces a quad face in the reversed orientation
  cells[1].vertices = {3u, 4u, 5u};
  cells[2].vertices = {1u, 3u, 2u, 5u};

  Triangulation<2> tria;

  tria.create_triangulation(vertices, cells, SubCellData());

  auto print_tria = [&]() {
    for (const auto &cell : tria.active_cell_iterators())
      {
        deallog << "cell = " << cell << std::endl;
        deallog << "  reference cell = " << cell->reference_cell().to_string()
                << std::endl;
        deallog << "  measure = " << cell->measure() << std::endl;
        for (unsigned int face_no : cell->face_indices())
          {
            deallog << "  face orientation " << face_no << " = "
                    << int(cell->combined_face_orientation(face_no))
                    << std::endl;
          }
        for (unsigned int vertex_no : cell->vertex_indices())
          {
            deallog << "  vertex " << vertex_no << " = "
                    << cell->vertex(vertex_no) << '\n';
          }
        for (unsigned int vertex_no : cell->vertex_indices())
          {
            deallog << "  vertex index " << vertex_no << " = "
                    << cell->vertex_index(vertex_no) << '\n';
          }
      }
  };

  print_tria();

  tria.refine_global(1);
  print_tria();

  tria.refine_global(1);
  print_tria();

#if 0
  std::ofstream out("out.vtu");

  GridOut().write_vtu(tria, out);
#endif
}
