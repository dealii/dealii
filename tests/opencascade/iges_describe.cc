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


// Read goteborg.iges and dump its topological structure to the
// logfile.

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

  TopoDS_Shape sh = read_IGES(SOURCE_DIR "/iges_files/goteborg.iges");
  std::vector<TopoDS_Face>   faces;
  std::vector<TopoDS_Edge>   edges;
  std::vector<TopoDS_Vertex> vertices;

  extract_geometrical_shapes(sh, faces, edges, vertices);

  std::tuple<unsigned int, unsigned int, unsigned int> n  = count_elements(sh);
  unsigned int                                         nf = std::get<0>(n);
  unsigned int                                         ne = std::get<1>(n);
  unsigned int                                         nv = std::get<2>(n);

  deallog << "Shape contains " << std::get<0>(n) << " faces, " << std::get<1>(n)
          << " edges, and " << std::get<2>(n) << " vertices." << std::endl;

  if (nf != faces.size())
    deallog << "Error!" << std::endl;
  if (ne != edges.size())
    deallog << "Error!" << std::endl;
  if (nv != vertices.size())
    deallog << "Error!" << std::endl;

  return 0;
}
