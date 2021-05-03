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


// Read a file in stl format, and use it to refine a grid. Then it
// verifies the final output

#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>

#include <deal.II/opencascade/manifold_lib.h>
#include <deal.II/opencascade/utilities.h>

#include <Standard_Stream.hxx>
#include <TopTools.hxx>
#include <TopoDS_Shape.hxx>

#include "../tests.h"

using namespace OpenCASCADE;

int
main()
{
  TopoDS_Shape        sh = read_STL(SOURCE_DIR "/stl_files/sphere_refined.stl");
  Triangulation<2, 3> tria;
  GridIn<2, 3>        gridin;
  gridin.attach_triangulation(tria);
  std::ifstream in;
  in.open(SOURCE_DIR "/stl_files/coarse_sphere_stl.inp");
  gridin.read_ucd(in, true);
  // OpenCASCADE::create_triangulation(shape, tria);
  OpenCASCADE::NormalToMeshProjectionManifold<2, 3> manifold(sh, 1e-3);
  // CachedNormalToMeshProjectionManifold manifold(shape,1e-3);
  tria.set_all_manifold_ids(0);
  tria.set_manifold(0, manifold);
  tria.refine_global(2);

  GridOut       gridout;
  std::ofstream ofile("output");
  gridout.write_ucd(tria, ofile);
}
