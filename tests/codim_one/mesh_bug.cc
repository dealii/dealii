
//----------------------------  mesh_bug.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2005 by the deal.II authors 
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  mesh_bug.cc  ---------------------------


// a short (a few lines) description of what the program does

#include "../tests.h"
#include <fstream>
#include <base/logstream.h>

// all include files you need here


#include <grid/tria.h>
#include <grid/tria_accessor.h>
#include <grid/grid_out.h>
#include <grid/grid_in.h>
#include <grid/tria_iterator.h>
#include <grid/grid_generator.h>
#include <grid/tria_boundary_lib.h>


int main () 
{
  std::ofstream logfile("mesh_bug/output");

  Triangulation<2,3> tria;
  GridIn<2,3> gi;
  gi.attach_triangulation(tria);
  std::ifstream infile("mesh_bug/cmp/generic");
  gi.read(infile);

  const std::vector<Point<3> > &vertices = tria.get_vertices();

  for(unsigned int i=0; i<vertices.size(); ++i)
    if(vertices[i](2)>1e-7)
      std::cout << "Error!" << std::endl;
  
  
  GridOut go;
  go.write_ucd(tria, logfile);

  return 0;
}
                  
