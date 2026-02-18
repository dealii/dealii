// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2014 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


// Read a file in iges format, and write it out again in the same
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
  TopoDS_Shape sh = read_IGES(SOURCE_DIR "/iges_files/wigley.iges");
  write_IGES(sh, "tmp.iges");
  std::ifstream in("tmp.iges");
  std::ofstream out("output");
  std::string   line;
  unsigned int  counter = 5;
  while (counter--)
    std::getline(in, line);
  while (std::getline(in, line))
    out << line << std::endl;
  in.close();
  out.close();
  return 0;
}
