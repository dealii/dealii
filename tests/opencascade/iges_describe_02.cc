// ---------------------------------------------------------------------
//
// Copyright (C) 2014 - 2018 by the deal.II authors
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
  std::vector<TopoDS_Compound>  compounds;
  std::vector<TopoDS_CompSolid> compsolids;
  std::vector<TopoDS_Solid>     solids;
  std::vector<TopoDS_Shell>     shells;
  std::vector<TopoDS_Wire>      wires;
  extract_compound_shapes(sh, compounds, compsolids, solids, shells, wires);

  deallog << "Shape contains " << compounds.size() << " compounds, "
          << compsolids.size() << " compound solids, " << solids.size()
          << " solids, " << shells.size() << " shells, " << wires.size()
          << " wires." << std::endl;
  return 0;
}
