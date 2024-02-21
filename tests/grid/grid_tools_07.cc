// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

//
// check find_active_cell_around_point finds only marked cells.
// If there are no marked cell containing the point it should not find
// anything. Latter is especially important if run with MPI.
//

#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_tools_cache.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"

/*
 * Generate a non-matching grid consisting of two elements which
 * do not share any face
 */

void
generate_grid(Triangulation<2> &triangulation)
{
  //   +-------------------+     +-------------------+
  //   |         |         |     |                   |
  //   |    2    |    3    |     |                   |
  //   |         |         |     |                   |
  // 10|_________|_________|8   9|        4          |
  //   |         |         |     |                   |
  //   |         |         |     |                   |
  //   |    0    |    1    |  p  |                   |
  //   |         |         |     |                   |
  //   +-------------------+     +-------------------+

  Triangulation<2> triangulation0;
  Triangulation<2> triangulation1;
  GridGenerator::hyper_rectangle(triangulation0, {-.5, -.5}, {.5, .5}, true);
  for (const auto &face : triangulation0.active_face_iterators())
    if (face->at_boundary() && face->boundary_id() == 1)
      face->set_boundary_id(8);
  for (const auto &face : triangulation0.active_face_iterators())
    if (face->at_boundary() && face->boundary_id() == 0)
      face->set_boundary_id(10);
  triangulation0.refine_global(1);
  GridGenerator::flatten_triangulation(triangulation0, triangulation1);

  Triangulation<2> triangulation2;
  GridGenerator::hyper_rectangle(triangulation2, {.5, -.5}, {1.5, .5}, true);
  for (const auto &face : triangulation2.active_face_iterators())
    if (face->at_boundary() && face->boundary_id() == 0)
      face->set_boundary_id(9);

  // tolerance 0. to ensure vertices are not merged
  GridGenerator::merge_triangulations(
    triangulation1, triangulation2, triangulation, 0., false, true);
}

std::vector<unsigned int>
mark_vertices_at_boundary(const types::boundary_id boundary_id,
                          const Triangulation<2>  &triangulation,
                          std::vector<bool>       &marked_vertices)
{
  std::vector<unsigned int> marked_cell_idxs;

  // mark_vertices_of_cells_at_boundary
  for (const auto &cell : triangulation.active_cell_iterators())
    {
      for (unsigned int face = 0; face < cell->n_faces(); ++face)
        {
          if (cell->face(face)->at_boundary() &&
              cell->face(face)->boundary_id() == boundary_id)
            {
              marked_cell_idxs.push_back(cell->index());
              for (unsigned int v = 0; v < cell->face(face)->n_vertices(); ++v)
                {
                  marked_vertices[cell->face(face)->vertex_index(v)] = true;
                }
              continue;
            }
        }
    }
  return marked_cell_idxs;
}

int
find_cell_at_point(const Point<2>          &p,
                   const Triangulation<2>  &triangulation,
                   const std::vector<bool> &marked_vertices,
                   const Triangulation<2>::active_cell_iterator &cell_hint,
                   const double tolerance = 1e-6)
{
  const MappingQ<2>            mapping{1};
  const GridTools::Cache<2, 2> cache{triangulation, mapping};
  const auto cell_and_pnt = GridTools::find_active_cell_around_point(
    cache.get_mapping(),
    cache.get_triangulation(),
    p,
    cache.get_vertex_to_cell_map(),
    cache.get_vertex_to_cell_centers_directions(),
    cell_hint,
    marked_vertices,
    cache.get_used_vertices_rtree(),
    tolerance,
    &cache.get_locally_owned_cell_bounding_boxes_rtree());

  if (cell_and_pnt.first != triangulation.active_cell_iterators().end())
    return cell_and_pnt.first->index();
  return -1;
}

int
main()
{
  // Setup
  Triangulation<2> triangulation;
  generate_grid(triangulation);

  const Point<2> p{.5, -.25}; // point located between cell 4 and 1
  const auto     cell_hint = ++triangulation.begin_active(); // hint cell 1

  initlog();
  deallog << std::setprecision(4);

  // Test1
  deallog << "Test if cell not found:\n\n";
  deallog << "Searching for point with marked cells that do not hold the "
             "point:\n";
  {
    std::vector<bool> marked_vertices(triangulation.n_vertices(), false);
    const auto        cell_idxs =
      mark_vertices_at_boundary(10, triangulation, marked_vertices);
    const int found_cell =
      find_cell_at_point(p, triangulation, marked_vertices, cell_hint);

    deallog << "Marked cells ";
    for (const auto &cell_idx : cell_idxs)
      deallog << cell_idx << " ";
    deallog << "\nFound cell: " << found_cell << "\n";
    deallog << "\n---------------\n\n";
  }

  // Test2
  deallog << "Test if marked cell is found with cell hint that has the "
             "point but is not marked:\n\n";
  {
    std::vector<bool> marked_vertices(triangulation.n_vertices(), false);

    // mark cell 4
    const auto cell_idxs =
      mark_vertices_at_boundary(9, triangulation, marked_vertices);
    deallog << "Searching for cells ";
    for (const auto &cell_idx : cell_idxs)
      deallog << cell_idx << " ";
    deallog << "with cell_hint " << cell_hint->index() << ":\n";

    const int found_cell =
      find_cell_at_point(p, triangulation, marked_vertices, cell_hint);
    deallog << "Found cell: " << found_cell << "\n";
  }
  deallog << std::endl << std::flush;

  return 0;
}
