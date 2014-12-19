// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


// similar to the hyper_ball_3d test, but for the cylinder grid. here,
// the cause for the failure was different, though: the description of
// the cells was wrong, and they were not sent through the
// GridReordering class which would have cured the problem.


#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/grid_generator.h>
#include <fstream>
#include <iomanip>



int main ()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);
  deallog << std::setprecision (2);

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
