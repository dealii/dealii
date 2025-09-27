// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2014 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Create a bspline, and write it as an IGES file.

#include <deal.II/opencascade/utilities.h>

#include <Standard_Stream.hxx>
#include <TopTools.hxx>
#include <TopoDS_Shape.hxx>

#include "../tests.h"

using namespace OpenCASCADE;

int
main()
{
  // Create a bspline passing through the points
  std::vector<Point<3>> pts;
  pts.push_back(Point<3>(0, 0, 0));
  pts.push_back(Point<3>(0, 1, 0));
  pts.push_back(Point<3>(1, 1, 0));
  pts.push_back(Point<3>(1, 0, 0));

  TopoDS_Edge edge = interpolation_curve(pts);
  write_IGES(edge, "tmp.iges");
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
