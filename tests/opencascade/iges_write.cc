//-----------------------------------------------------------
//
//    Copyright (C) 2014 by the deal.II authors 
//
//    This file is subject to LGPL and may not be distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-----------------------------------------------------------

// Read a file in iges format, and write it out again in the same
// format.

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
  TopoDS_Shape sh = read_IGES(SOURCE_DIR "/iges_files/wigley.iges");
  unsigned int nf=0, ne=0, nv=0;
  count_elements(sh, nf, ne, nv);

  write_IGES(sh, "output");
  return 0;
}
