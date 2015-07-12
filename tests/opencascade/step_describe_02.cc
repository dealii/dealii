//-----------------------------------------------------------
//
//    Copyright (C) 2014 - 2015 by the deal.II authors 
//
//    This file is subject to LGPL and may not be distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-----------------------------------------------------------

// Read goteborg.iges and dump its topological structure to the
// logfile.

#include "../tests.h"
#include <fstream>
#include <deal.II/base/logstream.h>

#include <deal.II/opencascade/utilities.h>
#include <TopTools.hxx>
#include <TopoDS_Shape.hxx>
#include <Standard_Stream.hxx>

using namespace OpenCASCADE;

int main () 
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  
  TopoDS_Shape sh = read_STEP(SOURCE_DIR "/step_files/goteborg.step");
  std::vector<TopoDS_Compound> compounds;
  std::vector<TopoDS_CompSolid> compsolids;
  std::vector<TopoDS_Solid> solids;
  std::vector<TopoDS_Shell> shells;
  std::vector<TopoDS_Wire> wires;
  extract_compound_shapes(sh, compounds, compsolids, solids, shells, wires);
  
  deallog << "Shape contains " << compounds.size() << " compounds, "
	  << compsolids.size() << " compound solids, "
	  << solids.size() << " solids, "
	  << shells.size() << " shells, "
	  << wires.size() << " wires." << std::endl;
  return 0;
}
                  
