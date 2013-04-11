//----------------------------  vertex_as_face_11.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2010, 2011, 2013 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  vertex_as_face_11.cc  ---------------------------

// verify that we can do things like cell->face() in 1d as well. here:
// check that we can default construct such an iterator and then
// assign to it using operator= (default construction did not work originally)


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

  typename Triangulation<1,spacedim>::face_iterator face;
  face = tria.begin_active()->face(0);

  deallog << face << std::endl;
}



int main ()
{
  std::ofstream logfile("vertex_as_face_13/output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  test<1> ();
  test<2> ();

  return 0;
}
