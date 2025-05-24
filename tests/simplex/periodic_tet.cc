/* ------------------------------------------------------------------------
 *
 * SPDX-License-Identifier: LGPL-2.1-or-later
 * Copyright (C) 2020 - 2024 by the deal.II authors
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

// Test tetrahedra + periodicity.

#include <deal.II/base/types.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/mapping_fe.h>

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_orientation.h>

#include <deal.II/lac/affine_constraints.h>

#include "../tests.h"

namespace
{
  void
  subdivided_hyper_rectangle_with_simplices_periodic(
    Triangulation<3>                &tria,
    const std::vector<unsigned int> &repetitions,
    const Point<3>                  &p1,
    const Point<3>                  &p2,
    const bool                       colorize)
  {
    std::vector<Point<3>>    vertices;
    std::vector<CellData<3>> cells;

    // determine cell sizes
    const Point<3> dx((p2[0] - p1[0]) / repetitions[0],
                      (p2[1] - p1[1]) / repetitions[1],
                      (p2[2] - p1[2]) / repetitions[2]);

    // create vertices
    for (unsigned int k = 0; k <= repetitions[2]; ++k)
      for (unsigned int j = 0; j <= repetitions[1]; ++j)
        for (unsigned int i = 0; i <= repetitions[0]; ++i)
          vertices.push_back(
            Point<3>(p1[0] + dx[0] * i, p1[1] + dx[1] * j, p1[2] + dx[2] * k));

    // create cells
    for (unsigned int k = 0; k < repetitions[2]; ++k)
      for (unsigned int j = 0; j < repetitions[1]; ++j)
        for (unsigned int i = 0; i < repetitions[0]; ++i)
          {
            // create reference HEX cell
            std::array<unsigned int, 8> quad{
              {(k + 0) * (repetitions[0] + 1) * (repetitions[1] + 1) +
                 (j + 0) * (repetitions[0] + 1) + i + 0,
               (k + 0) * (repetitions[0] + 1) * (repetitions[1] + 1) +
                 (j + 0) * (repetitions[0] + 1) + i + 1,
               (k + 0) * (repetitions[0] + 1) * (repetitions[1] + 1) +
                 (j + 1) * (repetitions[0] + 1) + i + 0,
               (k + 0) * (repetitions[0] + 1) * (repetitions[1] + 1) +
                 (j + 1) * (repetitions[0] + 1) + i + 1,
               (k + 1) * (repetitions[0] + 1) * (repetitions[1] + 1) +
                 (j + 0) * (repetitions[0] + 1) + i + 0,
               (k + 1) * (repetitions[0] + 1) * (repetitions[1] + 1) +
                 (j + 0) * (repetitions[0] + 1) + i + 1,
               (k + 1) * (repetitions[0] + 1) * (repetitions[1] + 1) +
                 (j + 1) * (repetitions[0] + 1) + i + 0,
               (k + 1) * (repetitions[0] + 1) * (repetitions[1] + 1) +
                 (j + 1) * (repetitions[0] + 1) + i + 1}};

            // TET cell 0
            {
              CellData<3> cell;
              cell.vertices = {{quad[0], quad[1], quad[3], quad[7]}};
              cells.push_back(cell);
            }

            // TET cell 1
            {
              CellData<3> cell;
              cell.vertices = {{quad[0], quad[1], quad[7], quad[5]}};
              cells.push_back(cell);
            }

            // TET cell 2
            {
              CellData<3> cell;
              cell.vertices = {{quad[0], quad[7], quad[3], quad[2]}};
              cells.push_back(cell);
            }

            // TET cell 3
            {
              CellData<3> cell;
              cell.vertices = {{quad[2], quad[6], quad[0], quad[7]}};
              cells.push_back(cell);
            }

            // TET cell 4
            {
              CellData<3> cell;
              cell.vertices = {{quad[4], quad[7], quad[5], quad[0]}};
              cells.push_back(cell);
            }

            // TET cell 5
            {
              CellData<3> cell;
              cell.vertices = {{quad[4], quad[6], quad[7], quad[0]}};
              cells.push_back(cell);
            }
          }

    // actually create triangulation
    tria.create_triangulation(vertices, cells, SubCellData());

    if (colorize)
      {
        // to colorize, run through all faces of all cells and set boundary
        // indicator to the correct value if it was 0.

        // use a large epsilon to compare numbers to avoid roundoff problems.
        double epsilon = std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < 3; ++i)
          epsilon = std::min(epsilon,
                             0.01 * (std::abs(p2[i] - p1[i]) / repetitions[i]));
        Assert(epsilon > 0,
               ExcMessage(
                 "The distance between corner points must be positive."));

        // run through all faces and check
        // if one of their center coordinates matches
        // one of the corner points. Comparisons
        // are made using an epsilon which
        // should be smaller than the smallest cell
        // diameter.
        auto face = tria.begin_face(), endface = tria.end_face();
        for (; face != endface; ++face)
          if (face->at_boundary())
            if (face->boundary_id() == 0)
              {
                const Point<3> center(face->center());

                if (std::abs(center[0] - p1[0]) < epsilon)
                  face->set_boundary_id(0);
                else if (std::abs(center[0] - p2[0]) < epsilon)
                  face->set_boundary_id(1);
                else if (std::abs(center[1] - p1[1]) < epsilon)
                  face->set_boundary_id(2);
                else if (std::abs(center[1] - p2[1]) < epsilon)
                  face->set_boundary_id(3);
                else if (std::abs(center[2] - p1[2]) < epsilon)
                  face->set_boundary_id(4);
                else if (std::abs(center[2] - p2[2]) < epsilon)
                  face->set_boundary_id(5);
                else
                  // triangulation says it
                  // is on the boundary,
                  // but we could not find
                  // on which boundary.
                  DEAL_II_ASSERT_UNREACHABLE();
              }
      }
  }
} // namespace

void
check(Triangulation<3> &tria)
{
  if (tria.n_levels() != 1)
    {
      deallog << "FAIL! Mesh has more than one refinement level." << std::endl;
    }

  // Check that each cell has a positive measure.
  {
    for (const auto &cell : tria.active_cell_iterators())
      {
        double measure = cell->measure();
        if (measure <= 0.0)
          {
            deallog << "FAIL! Found a cell with nonpositive measure: "
                    << cell->active_cell_index() << std::endl;
          }
      }
  }

  // Collect periodic face pairs and ensure that they make sense and we have
  // the right numbers of them.
  std::array<std::vector<GridTools::PeriodicFacePair<
               typename Triangulation<3>::cell_iterator>>,
             3>
    pairs;
  GridTools::collect_periodic_faces(tria, 0, 1, 0, pairs[0]);
  GridTools::collect_periodic_faces(tria, 2, 3, 1, pairs[1]);
  GridTools::collect_periodic_faces(tria, 4, 5, 2, pairs[2]);
  tria.add_periodicity(pairs[0]);
  tria.add_periodicity(pairs[1]);
  tria.add_periodicity(pairs[2]);

  {
    // Check that the vertices we found match in all but one coordinate
    // direction:
    unsigned int changing_direction = 0;
    for (const auto &pair : pairs)
      {
        for (const auto &face_pair : pair)
          {
            std::array<Point<3>, 3> face_1_vertices, face_2_vertices;
            for (unsigned int i = 0; i < 3; ++i)
              {
                face_1_vertices[i] =
                  face_pair.cell[0]->face(face_pair.face_idx[0])->vertex(i);
                face_2_vertices[i] =
                  face_pair.cell[1]->face(face_pair.face_idx[1])->vertex(i);
              }
            const auto permuted_face_1_vertices =
              ReferenceCells::Triangle.permute_by_combined_orientation(
                make_const_array_view(face_1_vertices), face_pair.orientation);

            deallog << "periodic face_pair:" << std::endl
                    << "    first  = (" << face_pair.cell[0]->index() << ", "
                    << face_pair.face_idx[0] << ")" << std::endl
                    << "    second = (" << face_pair.cell[1]->index() << ", "
                    << face_pair.face_idx[1] << ")" << std::endl
                    << "    relative orientation = "
                    << int(face_pair.orientation) << std::endl;

            bool fail = false;
            for (unsigned int i = 0; i < 3; ++i)
              for (unsigned int d = 0; d < 3; ++d)
                if (d != changing_direction)
                  fail = fail || std::abs(permuted_face_1_vertices[i][d] -
                                          face_2_vertices[i][d]) > 1e-10;

            if (fail)
              {
                deallog << "FAIL! Found non matching faces." << std::endl;
                return;
              }
          }
        ++changing_direction;
      }
  }

  // If we got this far, we have a mesh which supports periodicity.
  // Let's see if we can prescribe periodicity with actual DoFs now.
  // Check for C0 tetrahedra up to degree 3. (Highest implemented at the time of
  // writing.)
  for (const unsigned int degree : {1, 2, 3})
    {
      DoFHandler<3>  dof_handler(tria);
      FE_SimplexP<3> fe(degree);
      dof_handler.distribute_dofs(fe);
      std::array<AffineConstraints<double>, 3> constraints;

      // Periodic in the x direction:
      DoFTools::make_periodicity_constraints(
        dof_handler, 0, 1, 0, constraints[0]);
      constraints[0].close();

      // Periodic in the y direction:
      DoFTools::make_periodicity_constraints(
        dof_handler, 2, 3, 1, constraints[1]);
      constraints[1].close();

      // Periodic in the z direction:
      DoFTools::make_periodicity_constraints(
        dof_handler, 4, 5, 2, constraints[2]);
      constraints[2].close();

      // Some information for checking constraints.
      MappingFE<3> mapping(FE_SimplexP<3>(1));
      auto map = DoFTools::map_dofs_to_support_points(mapping, dof_handler);
      unsigned int changing_direction = 0;

      for (const auto &cons : constraints)
        {
          // Check the total number of constraints.
          deallog << "Number of constraints = " << cons.n_constraints()
                  << std::endl;

          // Loop through the constraints and verify that the correct components
          // of the DoF's support points are equal.
          for (const auto &constraint : cons.get_lines())
            {
              if (constraint.entries.size() == 1)
                {
                  if (std::abs(constraint.entries[0].second - 1.0) > 1e-10)
                    {
                      deallog << "FAIL! All constraint entries should be one."
                              << std::endl;
                      continue;
                    }

                  const types::global_dof_index dof_1 = constraint.index;
                  const types::global_dof_index dof_2 =
                    constraint.entries[0].first;
                  const std::array<bool, 3> results = {
                    {std::abs(map[dof_1][0] - map[dof_2][0]) > 1e-10,
                     std::abs(map[dof_1][1] - map[dof_2][1]) > 1e-10,
                     std::abs(map[dof_1][2] - map[dof_2][2]) > 1e-10}};
                  bool fail = false;

                  for (unsigned int i = 0; i < 3; ++i)
                    {
                      if (i == changing_direction)
                        {
                          continue;
                        }
                      else
                        {
                          fail = fail || results[i];
                        }
                    }

                  if (fail)
                    {
                      deallog << "FAIL! Found non matching constraint."
                              << std::endl;
                      deallog << "(" << map[dof_1][0] << "," << map[dof_1][1]
                              << "," << map[dof_1][2] << ")" << std::endl;
                      deallog << "(" << map[dof_2][0] << "," << map[dof_2][1]
                              << "," << map[dof_2][2] << ")" << std::endl;
                    }
                }
              else
                {
                  deallog
                    << "FAIL! "
                    << "Found constraint with the wrong number of entries."
                    << std::endl;
                  return;
                }
            }
          ++changing_direction;

          // If degree is 3, put down a little more data for the output.
          if (degree == 3)
            {
              cons.print(deallog.get_file_stream());
            }
        }
    }

  deallog << "OK!" << std::endl;
}

int
main()
{
  initlog();

  const Point<3>                  corner_p1(-0.5, 1.0, 0.0);
  const Point<3>                  corner_p2(0.5, 2.0, 1.2);
  const std::vector<unsigned int> sub_divisions{2, 2, 2};

  {
    Triangulation<3> tria;
    subdivided_hyper_rectangle_with_simplices_periodic(
      tria, sub_divisions, corner_p1, corner_p2, true);

    check(tria);
  }

  {
    Triangulation<3> hex_tria, tria;

    GridGenerator::subdivided_hyper_rectangle(
      hex_tria, sub_divisions, corner_p1, corner_p2, true);
    GridGenerator::convert_hypercube_to_simplex_mesh(hex_tria, tria);

    check(tria);
  }
}
