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

// Test that the function
// GridGenerator::subdivided_hyper_rectable_with_simplices() can give a mesh
// which is suitible for periodicity with tetrahedra.

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

void
check()
{
  Triangulation<3>          tria;
  Point<3>                  corner_p1     = {-0.5, 1.0, 0.0};
  Point<3>                  corner_p2     = {0.5, 2.0, 1.2};
  std::vector<unsigned int> sub_divisions = {2, 2, 2};
  GridGenerator::subdivided_hyper_rectangle_with_simplices(
    tria, sub_divisions, corner_p1, corner_p2, true, true);

  // Check that we have made a grid with the expected number of cells and only
  // one level of refinement.
  if (tria.n_cells() != 48)
    {
      deallog << "FAIL! Wrong number of cells." << std::endl;
    }

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

  // To enable periodicity the easy way, we need each boundary pair to be in
  // default orientation. Check that this holds.
  {
    IteratorFilters::AtBoundary        filter;
    const types::geometric_orientation default_orientation =
      numbers::default_geometric_orientation;
    const types::geometric_orientation default_mirror_orientation =
      internal::combined_face_orientation(false, false, false);
    for (const auto &pair : pairs)
      {
        for (const auto &match : pair)
          {
            if (match.orientation != default_orientation &&
                match.orientation != default_mirror_orientation)
              {
                deallog << "FAIL! Found a face pair which is not in default "
                        << "orientation." << std::endl;
                deallog << "Orientation is "
                        << static_cast<unsigned int>(match.orientation)
                        << " vs default "
                        << static_cast<unsigned int>(default_orientation)
                        << std::endl;
                deallog << "Cells " << match.cell[0]->active_cell_index()
                        << " and " << match.cell[1]->active_cell_index()
                        << std::endl;
                deallog << "Faces " << match.face_idx[0] << " and "
                        << match.face_idx[1] << std::endl;
                return;
              }
          }
      }
  }

  {
    // Check that we have gotten all the periodicity.
    bool         mismatched_faces   = false;
    unsigned int changing_direction = 0;
    for (const auto &pair : pairs)
      {
        if (pair.size() != 8)
          {
            deallog << "FAIL! There should be 8 matches on each face."
                    << std::endl;
            return;
          }
        for (int face = 0; face < static_cast<int>(pair.size()); ++face)
          {
            Point<3> p1(
              pair[face].cell[0]->face(pair[face].face_idx[0])->center());
            Point<3> p2(
              pair[face].cell[1]->face(pair[face].face_idx[1])->center());
            const std::array<bool, 3> results = {
              {std::abs(p1[0] - p2[0]) > 1e-10,
               std::abs(p1[1] - p2[1]) > 1e-10,
               std::abs(p1[2] - p2[2]) > 1e-10}};
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
                deallog << "FAIL! Found non matching faces." << std::endl;
                return;
              }
          }
        ++changing_direction;
      }
  }

  // write 2 outputs (total mesh and only surface mesh)
  const auto grid_out = [](const auto &tria, const bool surface_mesh_only) {
    GridOutFlags::Vtk flags;

    if (surface_mesh_only)
      {
        flags.output_cells         = false;
        flags.output_faces         = true;
        flags.output_edges         = false;
        flags.output_only_relevant = false;
      }

    GridOut grid_out;
    grid_out.set_flags(flags);

    grid_out.write_vtk(tria, deallog.get_file_stream());
  };

  grid_out(tria, false); // total mesh
  grid_out(tria, true);  // only surface mesh

  // If we got this far, we have a mesh which supports periodicity.
  // Let's see if we can proscribe periodicity with actual DoFs now.
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
          if (cons.n_constraints() != (2 * degree + 1) * (2 * degree + 1))
            {
              deallog << "FAIL! Wrong number of constraints!" << std::endl;
            }

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

  check();
}
