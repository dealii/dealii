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


#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include <deal.II/opencascade/boundary_lib.h>
#include <deal.II/opencascade/utilities.h>

#include <Standard_Stream.hxx>
#include <TopTools.hxx>
#include <TopoDS_Shape.hxx>

#include "../tests.h"

// Create a Triangulation, interpolate its boundary points to a smooth
// BSpline, and use that as an arlength Boundary Descriptor.

using namespace OpenCASCADE;

int
main()
{
  std::ofstream logfile("output");

  // Create a bspline passign through the points
  std::vector<Point<2>> pts;
  pts.push_back(Point<2>(0, 0));
  pts.push_back(Point<2>(0, 1));
  pts.push_back(Point<2>(1, 1));
  pts.push_back(Point<2>(1, 0));

  TopoDS_Edge edge = interpolation_curve(pts, Tensor<1, 2>(), true);
  ArclengthProjectionLineManifold<2, 2> manifold(edge);

  Triangulation<2> tria;
  GridGenerator::hyper_cube(tria);

  tria.set_all_manifold_ids_on_boundary(1);
  tria.begin_active()->face(2)->set_manifold_id(numbers::invalid_manifold_id);
  tria.set_manifold(1, manifold);

  tria.refine_global(3);
  GridOut gridout;
  gridout.write_msh(tria, logfile);

  return 0;
}
