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

// Makes a OpenCASCADE cicrular arc, and project a few points onto it.


#include "../tests.h"

#include <deal.II/opencascade/utilities.h>

#include <deal.II/base/logstream.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_boundary_lib.h>

#include <gp_Pnt.hxx>
#include <gp_Dir.hxx>
#include <gp_Ax2.hxx>
#include <GC_MakeCircle.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <TopoDS_Edge.hxx>

using namespace OpenCASCADE;


int main ()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  // A unit circle
  gp_Dir z_axis(0.,0.,1.);
  gp_Pnt center(0.,0.,0.);
  gp_Ax2 axis(center, z_axis);
  Standard_Real radius(1.);

  Handle(Geom_Curve) circle = GC_MakeCircle(axis, radius);
  TopoDS_Edge edge = BRepBuilderAPI_MakeEdge(circle);


  // Now get a few points and project
  // them on the circle
  std::vector<Point<3> > points;

  points.push_back(Point<3>(3,0,0));
  points.push_back(Point<3>(0,3,0));
  // This one is tricky... the point
  // is not on the xy plane. If we
  // put it on the axis (x=0, y=0),
  // there should not exist a
  // projection. Any value for z
  // should give the same result.
  points.push_back(Point<3>(.1,0,3));
  points.push_back(Point<3>(.1,0,4));


  double u, v;
  TopoDS_Shape sh;
  for (unsigned int i=0; i<points.size(); ++i)
    {
      std_cxx11::tuple<Point<3>, TopoDS_Shape, double, double> ref =
        project_point_and_pull_back(edge, points[i]);

      Point<3> pp = std_cxx11::get<0>(ref);
      sh = std_cxx11::get<1>(ref);
      u = std_cxx11::get<2>(ref);
      v = std_cxx11::get<3>(ref);

      deallog << "Origin: " << points[i]
              << ", on unit circle: " << pp
              << ", with local coordinates (u, v): (" << u
              << ", " << v << ")" << std::endl;
    }
  return 0;
}

