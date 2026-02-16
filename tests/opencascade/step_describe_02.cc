// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2014 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


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

  TopoDS_Shape sh = read_STEP(SOURCE_DIR "/step_files/goteborg.step");
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
