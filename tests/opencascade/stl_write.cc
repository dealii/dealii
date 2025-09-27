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


// Read a file in stl format, and write it out again in the same
// format.

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

  initlog();
  TopoDS_Shape sh = read_STL(SOURCE_DIR "/stl_files/sphere_refined.stl");
  write_STL(sh, "tmp.stl", 0.001, false, 1e-6, false, 0.001, false);
  cat_file("tmp.stl");
}
