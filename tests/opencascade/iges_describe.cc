// ---------------------------------------------------------------------
//
// Copyright (C) 2014 - 2020 by the deal.II authors
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
