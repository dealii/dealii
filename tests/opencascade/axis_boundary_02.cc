// ---------------------------------------------------------------------
//
// Copyright (C) 2014 - 2018 by the deal.II authors
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


// Test the class DirectionalProjectionBoundary

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include <deal.II/opencascade/boundary_lib.h>
#include <deal.II/opencascade/utilities.h>

#include <BRepFill.hxx>
#include <Standard_Stream.hxx>
#include <TopTools.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shape.hxx>

#include "../tests.h"

using namespace OpenCASCADE;

int
main()
{
  std::ofstream logfile("output");

  // Create a bspline passing through the points
  std::vector<Point<3>> pts1, pts2;
  pts1.push_back(Point<3>(0, 0, 0));
  pts1.push_back(Point<3>(1, 0, 0));

  pts2.push_back(Point<3>(0, 1, 0));
  pts2.push_back(Point<3>(.5, 1, 1));
  pts2.push_back(Point<3>(1, 1, 0));

  TopoDS_Edge edge1 = interpolation_curve(pts1);
  TopoDS_Edge edge2 = interpolation_curve(pts2);

  TopoDS_Face face = BRepFill::Face(edge1, edge2);

  DirectionalProjectionBoundary<2, 3> manifold(face, Point<3>(0, 0, 1));

  Triangulation<2, 3> tria;
  GridGenerator::hyper_cube(tria);

  tria.begin_active()->set_all_manifold_ids(1);
  tria.set_manifold(1, manifold);

  tria.refine_global(2);
  GridOut gridout;
  gridout.write_msh(tria, logfile);

  return 0;
}
