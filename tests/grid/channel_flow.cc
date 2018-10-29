// ---------------------------------------------------------------------
//
// Copyright (C) 2018 by the deal.II authors
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

// Test GridGenerator::channel_with_cylinder(). The output generated here is
// similar to that generated in extrude_copy_manifold.

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"

template <int dim>
void
test()
{
  Triangulation<dim> triangulation;
  GridGenerator::channel_with_cylinder(triangulation, 0.03, 4, 2.0, true);

  std::map<types::manifold_id, std::vector<Point<dim>>>
    manifold_to_face_centers;
  std::map<types::boundary_id, std::vector<Point<dim>>>
    boundary_to_face_centers;

  for (const auto &face : triangulation.active_face_iterators())
    {
      manifold_to_face_centers[face->manifold_id()].push_back(face->center());
      boundary_to_face_centers[face->boundary_id()].push_back(face->center());
    }

  deallog << "face centers:" << std::endl;
  for (std::pair<const types::manifold_id, std::vector<Point<dim>>> &key :
       manifold_to_face_centers)
    {
      if (key.first == numbers::flat_manifold_id)
        continue;
      deallog << "manifold id: " << key.first << std::endl;
      for (const Point<dim> &point : key.second)
        {
          deallog << point[0] << ", " << point[1];
          if (dim == 3)
            deallog << ", " << point[2];
          deallog << std::endl;
        }
    }
  for (std::pair<const types::boundary_id, std::vector<Point<dim>>> &key :
       boundary_to_face_centers)
    {
      if (key.first == numbers::internal_face_boundary_id)
        continue;
      deallog << "boundary id: " << key.first << std::endl;
      for (const Point<dim> &point : key.second)
        {
          deallog << point[0] << ", " << point[1];
          if (dim == 3)
            deallog << ", " << point[2];
          deallog << std::endl;
        }
    }

  // reenable if we want to look at pictures
  if (false)
    {
      std::ofstream out("out-" + std::to_string(dim) + ".vtu");
      GridOut       grid_out;
      grid_out.write_vtu(triangulation, out);
    }
}


int
main()
{
  initlog();

  deallog.push("2d");
  test<2>();
  deallog.pop();

  deallog.push("3d");
  test<3>();
  deallog.pop();
}
