
//----------------------------  iges_read.cc  ---------------------------
//    $Id: testsuite.html 13373 2006-07-13 13:12:08Z kanschat $
//    Version: $Name$ 
//
//    Copyright (C) 2005 by the deal.II authors 
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  iges_read.cc  ---------------------------


// Read goteborg.iges and dump its topological structure to the logfile.

#include "../tests.h"
#include <fstream>
#include <base/logstream.h>

#include <deal.II/grid/occ_utilities.h>
#include <TopTools.hxx>
#include <TopoDS_Shape.hxx>
#include <Standard_Stream.hxx>

using namespace OpenCASCADE;

int main () 
{
  // std::ofstream logfile("output");
  // deallog.attach(logfile);
  // deallog.depth_console(0);
  
  TopoDS_Shape sh = read_IGES(SOURCE_DIR "/iges_files/wigley.iges");
  unsigned int nf=0, ne=0, nv=0;
  count_elements(sh, nf, ne, nv);
  // deallog << "Shape contains " << nf << " faces, "
  // 	  << ne << " edges, and " 
  // 	  << nv << " vertices." << std::endl;

  write_IGES(sh, "output");

  // try 
  //   {
  //     TopoDS_Shape sh_fail = read_IGES("iges_files/nonexistent.iges");
  //   }
  // catch(...)
  //   {
  //     deallog << "Missing file correctly threw an exception." << std::endl;
  //   }
  
  return 0;
}
                  
