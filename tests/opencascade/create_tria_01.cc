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


// Read the file side.iges, attach it to a projector,
// create a single-cell Triangulation, and refine it with the new
// projector.

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

  TopoDS_Shape sh = read_IGES(SOURCE_DIR "/iges_files/wigley.iges", 0.001);
  std::vector<TopoDS_Face>   faces;
  std::vector<TopoDS_Edge>   edges;
  std::vector<TopoDS_Vertex> vertices;

  extract_geometrical_shapes(sh, faces, edges, vertices);

  // Create a boundary projector on the first face.
  NormalProjectionManifold<2, 3> boundary(faces[0]);

  // Create a Triangulation with a single cell
  Triangulation<2, 3> tria;
  create_triangulation(faces[0], tria);

  // Set the boundary
  tria.set_manifold(1, boundary);
  tria.begin()->set_all_manifold_ids(1);

  tria.refine_global(2);

  // You can open the generated file with gmsh
  GridOut gridout;
  gridout.write_msh(tria, logfile);

  return 0;
}
