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
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------


// Create a BSpline surface, and test axis projection.

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
  std::vector<Point<3>> pts;
  pts.push_back(Point<3>(0, 0, 0));
  pts.push_back(Point<3>(1, 0, 0));
  pts.push_back(Point<3>(1, 0, -1));
  pts.push_back(Point<3>(0, 0, -1));

  TopoDS_Edge edge1 = interpolation_curve(pts);
  for (unsigned int i = 0; i < pts.size(); ++i)
    pts[i] += Point<3>(0, 1, 0);
  TopoDS_Edge edge2 = interpolation_curve(pts);

  TopoDS_Face face = BRepFill::Face(edge1, edge2);


  NormalToMeshProjectionManifold<2, 3> manifold(face);

  Triangulation<2, 3> tria;
  GridGenerator::hyper_cube(tria);

  tria.begin_active()->set_all_manifold_ids(1);
  tria.set_manifold(1, manifold);

  tria.refine_global(2);
  GridOut gridout;
  gridout.write_msh(tria, logfile);

  return 0;
}
