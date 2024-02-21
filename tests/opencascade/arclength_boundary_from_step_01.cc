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


#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include <deal.II/opencascade/manifold_lib.h>
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

  TopoDS_Shape new_edge = read_STEP(SOURCE_DIR "/step_files/line_03.step", 1.0);

  std::vector<TopoDS_Face>   faces;
  std::vector<TopoDS_Edge>   edges;
  std::vector<TopoDS_Vertex> vertices;

  extract_geometrical_shapes(new_edge, faces, edges, vertices);
  AssertDimension(edges.size(), 1);

  ArclengthProjectionLineManifold<1, 3> manifold(edges[0]);

  Triangulation<1, 3> tria;
  GridGenerator::hyper_cube(tria);

  tria.begin_active()->set_all_manifold_ids(1);
  tria.set_manifold(1, manifold);

  tria.refine_global(4);
  GridOut gridout;
  gridout.write_msh(tria, logfile);

  return 0;
}
