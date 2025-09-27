// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Tests for project_real_point_to_unit_point_on_face in 2D and 3D on a cube,
// and also in 3D with a parallelepiped, checks against expected values.

#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <iostream>

#include "../tests.h"

void
dim2_grid()
{
  Triangulation<2> triangulation;

  const Point<2> p_ll(-1, -1); // lower left point
  const Point<2> p_ur(1, 1);   // upper right point
  GridGenerator::hyper_rectangle(triangulation, p_ll, p_ur, false);

  const Point<2> testp(.5, -.5); // test point

  MappingQ<2> mapping(1);

  deallog << "Check project for 2D cube from (-1,-1), to (1,1)." << std::endl;

  Triangulation<2>::active_cell_iterator cell = triangulation.begin_active(),
                                         endc = triangulation.end();
  deallog << "The test point has real coordinates: " << testp
          << ", and unit coordinates: "
          << mapping.transform_real_to_unit_cell(cell, testp) << std::endl;
  for (; cell != endc; ++cell)
    for (const unsigned int face : GeometryInfo<2>::face_indices())
      {
        deallog << "For face: " << face << ", the test point projects to: "
                << mapping.project_real_point_to_unit_point_on_face(cell,
                                                                    face,
                                                                    testp)
                << std::endl;
      }
}
void
dim3_grid()
{
  Triangulation<3> triangulation;

  const Point<3> p_ll(-1, -1, -1); // lower left point
  const Point<3> p_ur(1, 1, 1);    // upper right point
  GridGenerator::hyper_rectangle(triangulation, p_ll, p_ur, false);

  const Point<3> testp(.5, -.5, 0); // test point

  MappingQ<3> mapping(1);

  deallog << "Check project for 3D cube from (-1,-1,-1) to (1,1,1)."
          << std::endl;

  Triangulation<3>::active_cell_iterator cell = triangulation.begin_active(),
                                         endc = triangulation.end();
  deallog << "The test point has real coordinates: " << testp
          << ", and unit coordinates: "
          << mapping.transform_real_to_unit_cell(cell, testp) << std::endl;
  for (; cell != endc; ++cell)
    for (const unsigned int face : GeometryInfo<3>::face_indices())
      {
        deallog << "For face: " << face << ", the test point projects to: "
                << mapping.project_real_point_to_unit_point_on_face(cell,
                                                                    face,
                                                                    testp)
                << std::endl;
      }
}
void
dim3_parallelepiped_grid()
{
  Triangulation<3> triangulation;

  Point<3>(corners)[3];
  corners[0] = Point<3>(2, 0, 0);
  corners[1] = Point<3>(0, 2, 0);
  corners[2] = Point<3>(0, 1, 2);

  GridGenerator::parallelepiped(triangulation, corners);

  const Point<3> testp(1, 1, 1); // test point

  MappingQ<3> mapping(1);

  deallog
    << "Check project for 3D parallelepiped with vectors (2, 0, 0), (0, 2, 0), and (0, 1, 2)."
    << std::endl;
  Triangulation<3>::active_cell_iterator cell = triangulation.begin_active(),
                                         endc = triangulation.end();
  deallog << "The test point has real coordinates: " << testp
          << ", and unit coordinates: "
          << mapping.transform_real_to_unit_cell(cell, testp) << std::endl;

  for (; cell != endc; ++cell)
    for (const unsigned int face : GeometryInfo<3>::face_indices())
      {
        deallog << "For face: " << face << ", the test point projects to: "
                << mapping.project_real_point_to_unit_point_on_face(cell,
                                                                    face,
                                                                    testp)
                << std::endl;
      }
}
int
main()
{
  initlog();
  deallog.precision(3);

  dim2_grid();
  dim3_grid();
  dim3_parallelepiped_grid();

  return 0;
}
