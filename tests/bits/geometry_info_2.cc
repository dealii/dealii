//----------------------------  geometry_info_2.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2003, 2004, 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  geometry_info_2.cc  ---------------------------


// output most all static members of GeometryInfo

#include "../tests.h"
#include <base/logstream.h>
#include <base/geometry_info.h>

#include <fstream>
#include <cstdlib>

template <int dim>
void test ()
{
  deallog << GeometryInfo<dim>::children_per_cell << std::endl;
  deallog << GeometryInfo<dim>::faces_per_cell << std::endl;
  deallog << GeometryInfo<dim>::subfaces_per_face << std::endl;
  deallog << GeometryInfo<dim>::vertices_per_cell << std::endl;

  deallog << GeometryInfo<dim>::vertices_per_face << std::endl;
  deallog << GeometryInfo<dim>::lines_per_face << std::endl;
  deallog << GeometryInfo<dim>::quads_per_face << std::endl;
  
  deallog << GeometryInfo<dim>::lines_per_cell << std::endl;
  deallog << GeometryInfo<dim>::quads_per_cell << std::endl;
  deallog << GeometryInfo<dim>::hexes_per_cell << std::endl;
}


int main () 
{
  std::ofstream logfile("geometry_info_2.output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test<1> ();
  test<2> ();
  test<3> ();

  return 0;
}
