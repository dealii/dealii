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

// Create a bspline, and write it as an IGES file.

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
  // Create a bspline passing through the points
  std::vector<Point<3> > pts;
  pts.push_back(Point<3>(0,0,0));
  pts.push_back(Point<3>(0,1,0));
  pts.push_back(Point<3>(1,1,0));
  pts.push_back(Point<3>(1,0,0));

  TopoDS_Edge edge = interpolation_curve(pts);
  write_STEP(edge, "tmp.step");
  std::ifstream in("tmp.step");
  std::ofstream out("output");
  std::string line;
  unsigned int counter = 0;

  while (std::getline(in,line))
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

