// ---------------------------------------------------------------------
//
// Copyright (C) 2015 by the deal.II Authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

// Tests for project_real_point_to_unit_point_on_face in 2D and 3D on a cube,
// and also in 3D with a parallelepiped, checks against expected values.

#include "../tests.h"
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/fe/mapping_q1.h>
#include <iostream>
#include <fstream>
#include <cmath>

using namespace dealii;
void dim2_grid ()
{
  Triangulation<2> triangulation;

  const Point<2> p_ll (-1, -1);//lower left point
  const Point<2> p_ur (1, 1);//upper right point
  GridGenerator::hyper_rectangle (triangulation, p_ll, p_ur, false);

  const Point<2> testp (.5, -.5);  //test point

  MappingQGeneric<2> mapping(1);

  deallog<<"Check project for 2D cube from (-1,-1), to (1,1)."<<std::endl;

  Triangulation<2>::active_cell_iterator
  cell = triangulation.begin_active(),
  endc = triangulation.end();
  deallog<<"The test point has real coordinates: "<<testp
         <<", and unit coordinates: "
         <<mapping.transform_real_to_unit_cell(cell, testp)<<std::endl;
  for (; cell!=endc; ++cell)
    for (unsigned int face=0; face<GeometryInfo<2>::faces_per_cell; ++face)
      {
        deallog<<"For face: "<<face<<", the test point projects to: "
               <<mapping.project_real_point_to_unit_point_on_face(cell, face, testp)<<std::endl;
      }
}
void dim3_grid ()
{
  Triangulation<3> triangulation;

  const Point<3> p_ll (-1, -1, -1);//lower left point
  const Point<3> p_ur (1, 1, 1);//upper right point
  GridGenerator::hyper_rectangle (triangulation, p_ll, p_ur, false);

  const Point<3> testp (.5, -.5, 0);  //test point

  MappingQGeneric<3> mapping(1);

  deallog<<"Check project for 3D cube from (-1,-1,-1) to (1,1,1)."<<std::endl;

  Triangulation<3>::active_cell_iterator
  cell = triangulation.begin_active(),
  endc = triangulation.end();
  deallog<<"The test point has real coordinates: "<<testp
         <<", and unit coordinates: "
         <<mapping.transform_real_to_unit_cell(cell, testp)<<std::endl;
  for (; cell!=endc; ++cell)
    for (unsigned int face=0; face<GeometryInfo<3>::faces_per_cell; ++face)
      {
        deallog<<"For face: "<<face<<", the test point projects to: "
               <<mapping.project_real_point_to_unit_point_on_face(cell, face, testp)<<std::endl;
      }
}
void dim3_parallelepiped_grid ()
{
  Triangulation<3> triangulation;

  Point<3> (corners) [3];
  corners[0] = Point<3> (2, 0, 0);
  corners[1] = Point<3> (0, 2, 0);
  corners[2] = Point<3> (0, 1, 2);

  GridGenerator::parallelepiped (triangulation, corners);

  const Point<3> testp (1, 1, 1);  //test point

  MappingQGeneric<3> mapping(1);

  deallog<<"Check project for 3D parallelepiped with vectors (2, 0, 0), (0, 2, 0), and (0, 1, 2)."<<std::endl;
  Triangulation<3>::active_cell_iterator
  cell = triangulation.begin_active(),
  endc = triangulation.end();
  deallog<<"The test point has real coordinates: "<<testp
         <<", and unit coordinates: "
         <<mapping.transform_real_to_unit_cell(cell, testp)<<std::endl;

  for (; cell!=endc; ++cell)
    for (unsigned int face=0; face<GeometryInfo<3>::faces_per_cell; ++face)
      {
        deallog<<"For face: "<<face<<", the test point projects to: "
               <<mapping.project_real_point_to_unit_point_on_face(cell, face, testp)<<std::endl;
      }
}
int main ()
{

  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.precision(3);

  dim2_grid ();
  dim3_grid ();
  dim3_parallelepiped_grid ();

  return 0;
}
