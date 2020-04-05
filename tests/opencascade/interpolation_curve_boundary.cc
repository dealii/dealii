// ---------------------------------------------------------------------
//
// Copyright (C) 2014 - 2019 by the deal.II authors
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


// Create a Triangulation, interpolate its boundary points to a close
// smooth BSpline, and use that as a Boundary Descriptor.

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/opencascade/boundary_lib.h>

#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <GC_MakeCircle.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Wire.hxx>
#include <gp_Ax2.hxx>
#include <gp_Dir.hxx>
#include <gp_Pnt.hxx>

#include "../tests.h"

using namespace OpenCASCADE;

int
main()
{
  initlog();

  // The curve passing through the vertices of the unit square.
  std::vector<Point<3>> vertices;
  vertices.push_back(Point<3>(0, 0, 0));
  vertices.push_back(Point<3>(1, 0, 0));
  vertices.push_back(Point<3>(1, 1, 0));
  vertices.push_back(Point<3>(0, 1, 0));
  TopoDS_Shape shape = interpolation_curve(vertices, Point<3>(), true);

  // Create a boundary projector.
  NormalProjectionBoundary<2, 3> boundary_line(shape);

  // The unit square.
  Triangulation<2, 3> tria;
  GridGenerator::hyper_cube(tria);

  // Set the exterior boundary
  tria.set_all_manifold_ids(0);
  tria.set_manifold(0, boundary_line);

  // This is here to ignore the points created in the interior of the
  // face.
  tria.begin()->set_manifold_id(1);

  // We refine twice, and expect the outer points to end up on a
  // smooth curve interpolating the square vertices.
  tria.refine_global(2);


  // You can open the generated file with gmsh.
  GridOut gridout;
  gridout.write_msh(tria, deallog.get_file_stream());

  return 0;
}
