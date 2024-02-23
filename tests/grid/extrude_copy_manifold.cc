// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Check that we can correctly copy manifold ids during extrusion. The output
// was manually checked by plotting points in paraview.


#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"

void
test()
{
  Triangulation<2> triangulation_2;
  GridGenerator::hyper_ball(triangulation_2);
  for (auto &cell : triangulation_2.active_cell_iterators())
    if (cell->at_boundary())
      {
        cell->set_manifold_id(1);
        for (const unsigned int face_n : GeometryInfo<2>::face_indices())
          if (!cell->face(face_n)->at_boundary())
            cell->face(face_n)->set_manifold_id(1);
      }
  for (auto &cell : triangulation_2.active_cell_iterators())
    if (!cell->at_boundary())
      cell->set_all_manifold_ids(numbers::flat_manifold_id);
  TransfiniteInterpolationManifold<2> tfi_manifold_2;
  tfi_manifold_2.initialize(triangulation_2);
  triangulation_2.set_manifold(1, tfi_manifold_2);

  Triangulation<3> triangulation_3;
  GridGenerator::extrude_triangulation(
    triangulation_2, 3, 1.0, triangulation_3, true);
  triangulation_3.set_manifold(0, FlatManifold<3>());
  triangulation_3.set_manifold(1, FlatManifold<3>());
  TransfiniteInterpolationManifold<3> tfi_manifold;
  tfi_manifold.initialize(triangulation_3);
  CylindricalManifold<3> cylinder_manifold(2);
  triangulation_3.set_manifold(0, cylinder_manifold);
  triangulation_3.set_manifold(42, tfi_manifold);
  triangulation_3.refine_global(1);

  std::map<types::manifold_id, std::vector<Point<3>>> line_centers;
  std::map<types::manifold_id, std::vector<Point<3>>> face_centers;
  std::map<types::manifold_id, std::vector<Point<3>>> cell_centers;

  for (const auto &cell : triangulation_3.active_cell_iterators())
    {
      cell_centers[cell->manifold_id()].push_back(cell->center());
      for (const unsigned int face_n : GeometryInfo<3>::face_indices())
        face_centers[cell->face(face_n)->manifold_id()].push_back(
          cell->face(face_n)->center());
      for (unsigned int line_n = 0; line_n < GeometryInfo<3>::lines_per_cell;
           ++line_n)
        line_centers[cell->line(line_n)->manifold_id()].push_back(
          cell->line(line_n)->center());
    }

  deallog << "cell centers:" << std::endl;
  for (const std::pair<const types::manifold_id, std::vector<Point<3>>> &key :
       cell_centers)
    {
      // cell centers are already unique: we only visited each cell once
      deallog << "manifold id: " << key.first << std::endl;
      for (const Point<3> &point : key.second)
        deallog << point[0] << ", " << point[1] << ", " << point[2]
                << std::endl;
    }
  deallog << std::endl;

  // The generated arrays of face and line centers contain duplicates: get rid
  // of them by sorting and then std::unique-ing
  auto point_comparator = [](const Point<3> &a, const Point<3> &b) {
    // to minimize roundoff problems, align numbers in some lexicographic-like
    // order
    return 1e-10 * a[0] + 1e-5 * a[1] + a[2] <
           1e-10 * b[0] + 1e-5 * b[1] + b[2];
  };

  deallog << "face centers:" << std::endl;
  for (std::pair<const types::manifold_id, std::vector<Point<3>>> &key :
       face_centers)
    {
      deallog << "manifold id: " << key.first << std::endl;
      std::sort(key.second.begin(), key.second.end(), point_comparator);
      key.second.erase(std::unique(key.second.begin(), key.second.end()),
                       key.second.end());
      for (const Point<3> &point : key.second)
        deallog << point[0] << ", " << point[1] << ", " << point[2]
                << std::endl;
    }

  deallog << "line centers:" << std::endl;
  for (std::pair<const types::manifold_id, std::vector<Point<3>>> &key :
       line_centers)
    {
      deallog << "manifold id: " << key.first << std::endl;
      std::sort(key.second.begin(), key.second.end(), point_comparator);
      key.second.erase(std::unique(key.second.begin(), key.second.end()),
                       key.second.end());
      for (const Point<3> &point : key.second)
        deallog << point[0] << ", " << point[1] << ", " << point[2]
                << std::endl;
    }

  // reenable if we want to look at pictures
  if (false)
    {
      std::ofstream out("out-" + std::to_string(3) + ".vtu");
      GridOut       grid_out;
      grid_out.write_vtu(triangulation_3, out);

      Triangulation<2, 3> triangulation_23;
      GridGenerator::extract_boundary_mesh(triangulation_3, triangulation_23);
      std::ofstream out23("out-23.vtu");
      grid_out.write_vtu(triangulation_23, out23);
    }
}


int
main()
{
  initlog();
  test();
}
