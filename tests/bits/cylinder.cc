//----------------------------  cylinder.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2003, 2004 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  cylinder.cc  ---------------------------

// similar to the hyper_ball_3d test, but for the cylinder grid. here,
// the cause for the failure was different, though: the description of
// the cells was wrong, and they were not sent through the
// GridReordering class which would have cured the problem.


#include "../tests.h"
#include <base/logstream.h>
#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <grid/tria_accessor.h>
#include <grid/grid_generator.h>
#include <fstream>

    

int main () 
{
  std::ofstream logfile("cylinder.output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  logfile.precision (2);

                                   // generate a cylinder
  Triangulation<3> tria;
  GridGenerator::cylinder (tria, std::sqrt(2.));

                                   // output all faces. here, we
                                   // should have 18 (two layers of
                                   // cells, each with 4 outer faces,
                                   // plus 5 faces each for the top
                                   // and bottom of the cylinder)
  unsigned int external_faces = 0;
  for (Triangulation<3>::face_iterator face=tria.begin_face();
       face!=tria.end_face(); ++face)
    {
      deallog << face << "   "
              << (int)face->boundary_indicator() << "  "
              << '<' << face->vertex(0) << '>' << std::endl
              << "           <" << face->vertex(1) << '>' << std::endl
              << "           <" << face->vertex(2) << '>' << std::endl
              << "           <" << face->vertex(3) << '>' << std::endl;
      if (face->at_boundary())
        ++external_faces;
    }

  deallog << "External faces: " << external_faces << std::endl;

  Assert (external_faces == 18, ExcInternalError());
  
  return 0;
}
