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
  TopoDS_Shape sh = read_STEP(SOURCE_DIR "/step_files/wigley.step", 1.0);
  write_STEP(sh, "tmp.step");
  std::ifstream in("tmp.step");
  std::ofstream out("output");
  std::string   line;
  unsigned int  counter = 0;

  while (std::getline(in, line))
    {
      counter++;
      if ((counter == 4) || (counter == 5) || (counter == 6) ||
          (counter == 18) || (counter == 19))
        {
        }
      else
        out << line << std::endl;
    }


  in.close();
  out.close();
  return 0;
}
