// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Create a NURBS with a quarter of a disc, make a Manifold out of it, create a
// Triangulation with it, attach the manifold, and make a few refinements

#include <deal.II/grid/grid_out.h>

#include <deal.II/opencascade/manifold_lib.h>
#include <deal.II/opencascade/utilities.h>

#include <BRepPrimAPI_MakeRevol.hxx>

#include "../tests.h"


using namespace OpenCASCADE;

int
main()
{
  std::ofstream logfile("output");

  // Create an edge, and make a face out of it by revolving along the z
  // axis.
  std::vector<Point<2>> ps;
  ps.push_back(Point<2>(0, 1));
  ps.push_back(Point<2>(0, 2));
  auto edge = interpolation_curve(ps);

  BRepPrimAPI_MakeRevol revol(edge,
                              gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1)),
                              -numbers::PI / 2);
  revol.Build();
  auto sh = revol.Shape();

  std::vector<TopoDS_Face>   faces;
  std::vector<TopoDS_Edge>   edges;
  std::vector<TopoDS_Vertex> vertices;

  extract_geometrical_shapes(sh, faces, edges, vertices);

  // Create a boundary projector on the first face.
  NURBSPatchManifold<2, 2> manifold(faces[0]);

  // Create a Triangulation with a single cell
  Triangulation<2, 2> tria;
  create_triangulation(faces[0], tria);

  // Set the manifold
  tria.set_manifold(1, manifold);
  tria.set_all_manifold_ids(1);

  tria.refine_global(3);

  // You can open the generated file with gmsh
  GridOut gridout;
  gridout.write_msh(tria, logfile);
}
