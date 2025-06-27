// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2014 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


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
  // This test might trigger spurious floating point exception despite
  // functioning properly. Simply disable floating point exceptions again
  // (after they had been enabled int tests.h)
#if defined(DEBUG) && defined(DEAL_II_HAVE_FP_EXCEPTIONS)
  {
    const int current_fe_except = fegetexcept();
    fedisableexcept(current_fe_except);
  }
#endif

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
