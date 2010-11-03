
//----------------------------  template.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2005, 2008, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  template.cc  ---------------------------


// read a couple meshes with a codimension and write them out again

#include "../tests.h"
#include <fstream>
#include <base/logstream.h>
#include <grid/tria.h>
#include <grid/grid_in.h>
#include <grid/grid_out.h>

#include <string>

std::ofstream logfile("grid_in_out/output");

template <int dim, int spacedim>
void test(std::string filename) {
    Triangulation<dim, spacedim> tria;
    GridIn<dim, spacedim> gi;
    gi.attach_triangulation (tria);
    std::ifstream in (filename.c_str());
    gi.read_ucd (in);

    GridOut grid_out;
    grid_out.set_flags (GridOutFlags::Ucd(true));
    grid_out.write_ucd (tria, logfile);
    grid_out.write_msh (tria, logfile);
}

int main ()
{
  deallog.attach(logfile);
  deallog.depth_console(0);

  test<1,2>("grids/circle_1.inp");
  test<2,3>("grids/square.inp");
  test<2,3>("grids/sphere_1.inp");

  return 0;
}

