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

// Read a file in iges format, and write it out again in the same
// format.

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
  TopoDS_Shape sh = read_STEP(SOURCE_DIR "/step_files/wigley.step");
  write_STEP(sh, "tmp.step");
  std::ifstream in("tmp.step");
  std::ofstream out("output");
  std::string line;
  unsigned int counter = 0;

  while(std::getline(in,line))
       {
       counter++;
       if ( (counter == 4) ||
            (counter == 5) ||
            (counter == 6) ||
            (counter == 18) ||
            (counter == 19)   )
          {
	  }
       else
          out << line << std::endl;
       }  
    
    
  in.close();
  out.close();
  return 0;
}
