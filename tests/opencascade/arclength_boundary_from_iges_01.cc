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

#include "../tests.h"
#include <fstream>
#include <deal.II/base/logstream.h>

#include <deal.II/opencascade/utilities.h>
#include <deal.II/opencascade/boundary_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_generator.h>

#include <TopTools.hxx>
#include <TopoDS_Shape.hxx>
#include <Standard_Stream.hxx>

// Create a Triangulation, interpolate its boundary points to a smooth
// BSpline, and use that as an arlength Boundary Descriptor.

using namespace OpenCASCADE;

int main ()
{
  std::ofstream logfile("output");

  TopoDS_Shape new_edge = read_IGES(SOURCE_DIR"/iges_files/line_03.iges", 1.0);

  ArclengthProjectionLineManifold<1,3> manifold(new_edge);

  Triangulation<1,3> tria;
  GridGenerator::hyper_cube(tria);

  tria.begin_active()->set_all_manifold_ids(1);
  tria.set_manifold(1, manifold);

  tria.refine_global(4);
  GridOut gridout;
  gridout.write_msh(tria, logfile);

  return 0;
}

