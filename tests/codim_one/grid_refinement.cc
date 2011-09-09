//----------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2005, 2008, 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------------------------------------------------


// see what happens when creating a surface mesh and then refining it

#include "../tests.h"
#include <fstream>
#include <deal.II/base/logstream.h>

// all include files you need here

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria_boundary_lib.h>

#include <fstream>
#include <string>

std::ofstream logfile("grid_refinement/output");

template <int dim, int spacedim>
void test(std::string filename) {
  HyperBallBoundary<dim, spacedim> boundary;
  Triangulation<dim, spacedim> tria;
  tria.set_boundary(1, boundary);
  GridIn<dim, spacedim> gi;
  gi.attach_triangulation (tria);
  std::ifstream in (filename.c_str());
  gi.read_ucd (in);

  GridOut grid_out;
  grid_out.set_flags (GridOutFlags::Ucd(true));
  for(unsigned int cycle=0; cycle<3; ++cycle) {
    tria.refine_global(1);
    grid_out.write_msh (tria, logfile);
  }
}

int main ()
{
  deallog.attach(logfile);
  deallog.depth_console(0);

  deallog << "Test<1,2>" << std::endl;
  test<1,2>("grids/circle_1.inp");

  deallog << std::endl << "Test<1,2>" << std::endl;
  test<2,3>("grids/sphere_1.inp");

  return 0;
}

