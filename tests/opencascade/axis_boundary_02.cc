// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2014 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test the class DirectionalProjectionManifold

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include <deal.II/opencascade/manifold_lib.h>
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

  DirectionalProjectionManifold<2, 3> manifold(face, Point<3>(0, 0, 1));

  Triangulation<2, 3> tria;
  GridGenerator::hyper_cube(tria);

  tria.begin_active()->set_all_manifold_ids(1);
  tria.set_manifold(1, manifold);

  tria.refine_global(2);
  GridOut gridout;
  gridout.write_msh(tria, logfile);

  return 0;
}
