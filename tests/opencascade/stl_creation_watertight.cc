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


// Create an STL starting from a watertight IGES file.

#include <deal.II/opencascade/utilities.h>

#include <Standard_Stream.hxx>
#include <TopTools.hxx>
#include <TopoDS_Shape.hxx>

#include "../tests.h"

using namespace OpenCASCADE;

int
main()
{
  initlog();

  const std::string closed_iges       = SOURCE_DIR "/iges_files/wing.iges";
  const std::string output_closed_stl = "output.stl";
  auto              shape_fuse = OpenCASCADE::read_IGES(closed_iges, 1e-2);
  write_STL(shape_fuse, output_closed_stl, 0.1, true, 1e-6, false, 0.1, false);
  deallog << "OK" << std::endl;
}
