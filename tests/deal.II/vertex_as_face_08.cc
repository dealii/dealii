//----------------------------  vertex_as_face_04.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2010, 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  vertex_as_face_04.cc  ---------------------------

// verify that we can do things like cell->face() in 1d as well. here:
// Triangulation::get_boundary_indicators()


#include "../tests.h"
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/grid_generator.h>

#include <fstream>


template <int spacedim>
void test ()
{
  Triangulation<1,spacedim> tria;
  GridGenerator::hyper_cube (tria);
  tria.refine_global (2);

  std::vector<unsigned char>
    boundary_ids = tria.get_boundary_indicators ();

  for (unsigned int i=0; i<boundary_ids.size(); ++i)
    deallog << (int)boundary_ids[i] << std::endl;
}



int main ()
{
  std::ofstream logfile("vertex_as_face_08/output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  test<1> ();
  test<2> ();

  return 0;
}
