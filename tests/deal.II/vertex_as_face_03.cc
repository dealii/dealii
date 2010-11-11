//----------------------------  vertex_as_face_03.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  vertex_as_face_03.cc  ---------------------------

// verify that we can do things like cell->face() in 1d as well. here:
// test cell->face()->at_boundary()


#include "../tests.h"
#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <grid/tria_accessor.h>
#include <grid/grid_generator.h>

#include <fstream>


template <int spacedim>
void test ()
{
  Triangulation<1,spacedim> tria;
  GridGenerator::hyper_cube (tria);

  deallog << "Coarse mesh:" << std::endl;
  deallog << "Left vertex=" << tria.begin_active()->face(0)->at_boundary() << std::endl;
  deallog << "Right vertex=" << tria.begin_active()->face(1)->at_boundary() << std::endl;

  tria.refine_global (2);

  for (typename Triangulation<1,spacedim>::active_cell_iterator
	 cell = tria.begin_active();
       cell != tria.end(); ++cell)
    {
      deallog << "Cell: " << cell << std::endl;
      deallog << "Left vertex=" << cell->face(0)->at_boundary() << std::endl;
      deallog << "Right vertex=" << cell->face(1)->at_boundary() << std::endl;
    }
}



int main ()
{
  std::ofstream logfile("vertex_as_face_03/output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  test<1> ();
  test<2> ();

  return 0;
}
