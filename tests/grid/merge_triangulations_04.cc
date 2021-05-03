// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2020 by the deal.II authors
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

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"


int
main()
{
  initlog();

  // set up the bulk triangulation
  Triangulation<2> bulk_triangulation;
  GridGenerator::subdivided_hyper_rectangle(bulk_triangulation,
                                            {22u, 4u},
                                            Point<2>(0.0, 0.0),
                                            Point<2>(2.2, 0.41));
  std::set<Triangulation<2>::active_cell_iterator> cells_to_remove;
  Tensor<1, 2> cylinder_triangulation_offset;
  for (const auto &cell : bulk_triangulation.active_cell_iterators())
    {
      if ((cell->center() - Point<2>(0.2, 0.2)).norm() < 0.15)
        cells_to_remove.insert(cell);

      if (cylinder_triangulation_offset == Point<2>())
        {
          for (const unsigned int vertex_n : GeometryInfo<2>::vertex_indices())
            if (cell->vertex(vertex_n) == Point<2>())
              {
                // skip two cells in the bottom left corner
                cylinder_triangulation_offset =
                  2.0 * (cell->vertex(3) - Point<2>());
                break;
              }
        }
    }
  Triangulation<2> result_1;
  GridGenerator::create_triangulation_with_removed_cells(bulk_triangulation,
                                                         cells_to_remove,
                                                         result_1);

  // set up the cylinder triangulation
  Triangulation<2> cylinder_triangulation;
  GridGenerator::hyper_cube_with_cylindrical_hole(cylinder_triangulation,
                                                  0.05,
                                                  0.41 / 4.0);
  GridTools::shift(cylinder_triangulation_offset, cylinder_triangulation);
  // dumb hack
  for (const auto &cell : cylinder_triangulation.active_cell_iterators())
    cell->set_material_id(2);

  // merge them together
  auto minimal_line_length = [](const Triangulation<2> &tria) -> double {
    double min_line_length = std::numeric_limits<double>::max();
    for (const auto &cell : tria.active_cell_iterators())
      for (unsigned int line_n = 0; line_n < GeometryInfo<2>::lines_per_cell;
           ++line_n)
        min_line_length = std::min(min_line_length,
                                   (cell->line(line_n)->vertex(0) -
                                    cell->line(line_n)->vertex(1))
                                     .norm());
    return min_line_length;
  };

  // the cylindrical triangulation might not match the Cartesian grid: as a
  // result the vertices might not be lined up. Get around this by deleting
  // duplicated vertices with a very low numerical tolerance.
  const double tolerance =
    std::min(minimal_line_length(result_1),
             minimal_line_length(cylinder_triangulation)) /
    2.0;


  Triangulation<2> result_2;
  GridGenerator::merge_triangulations(result_1,
                                      cylinder_triangulation,
                                      result_2,
                                      tolerance);

  const types::manifold_id tfi_id   = 1;
  const types::manifold_id polar_id = 0;
  for (const auto &cell : result_2.active_cell_iterators())
    {
      // set all non-boundary manifold ids to the new TFI manifold id.
      if (cell->material_id() == 2)
        {
          cell->set_manifold_id(tfi_id);
          for (const unsigned int face_n : GeometryInfo<2>::face_indices())
            {
              if (cell->face(face_n)->at_boundary())
                cell->face(face_n)->set_manifold_id(polar_id);
              else
                cell->face(face_n)->set_manifold_id(tfi_id);
            }
        }
    }

  PolarManifold<2> polar_manifold(Point<2>(0.2, 0.2));
  result_2.set_manifold(polar_id, polar_manifold);
  TransfiniteInterpolationManifold<2> inner_manifold;
  inner_manifold.initialize(result_2);
  result_2.set_manifold(tfi_id, inner_manifold);

  std::vector<Point<2> *> inner_pointers;
  for (const auto &cell : result_2.active_cell_iterators())
    {
      for (const unsigned int face_n : GeometryInfo<2>::face_indices())
        {
          if (cell->face(face_n)->manifold_id() == polar_id)
            {
              inner_pointers.push_back(&cell->face(face_n)->vertex(0));
              inner_pointers.push_back(&cell->face(face_n)->vertex(1));
            }
        }
    }
  // de-duplicate
  std::sort(inner_pointers.begin(), inner_pointers.end());
  inner_pointers.erase(std::unique(inner_pointers.begin(),
                                   inner_pointers.end()),
                       inner_pointers.end());

  // find the current center...
  Point<2> center;
  for (const Point<2> *const ptr : inner_pointers)
    center += *ptr / double(inner_pointers.size());

  // and recenter at (0.2, 0.2)
  for (Point<2> *const ptr : inner_pointers)
    *ptr += Point<2>(0.2, 0.2) - center;

  result_2.refine_global(1);
  for (const auto &cell : result_2.active_cell_iterators())
    {
      if (cell->at_boundary())
        {
          for (const unsigned int face_n : GeometryInfo<2>::face_indices())
            {
              auto face = cell->face(face_n);
              if (face->at_boundary())
                {
                  deallog << "boundary face with center " << face->center()
                          << std::endl;
                }
            }
        }
    }
  deallog << "OK" << std::endl;
}
