// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2014 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Create a circle, a Triangulation, and try to project normally on
// it.

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/opencascade/manifold_lib.h>

#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepPrimAPI_MakeSphere.hxx>
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
  std::ofstream logfile("output");

  gp_Pnt        center(.5, .5, .5);
  Standard_Real radius(Point<3>().distance(point<3>(center)));


  TopoDS_Face face = BRepPrimAPI_MakeSphere(center, radius);

  // Create a boundary projector.
  NormalProjectionManifold<3, 3> sphere(face);


  // The unit cube.
  Triangulation<3, 3> tria;
  GridGenerator::hyper_cube(tria);

  // Set the exterior boundary
  tria.set_manifold(1, sphere);

  // This is here to ignore the
  // points created in the interior
  // of the face.
  tria.begin()->set_all_manifold_ids(1);
  tria.begin()->set_manifold_id(0);
  tria.set_manifold(0, FlatManifold<3>());

  tria.refine_global(2);


  // You can open the generated file with gmsh
  GridOut gridout;
  gridout.write_msh(tria, logfile);

  return 0;
}
