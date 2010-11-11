//----------------------------  vertex_as_face_04.cc  ---------------------------
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
//----------------------------  vertex_as_face_04.cc  ---------------------------

// verify that we can do things like cell->face() in 1d as well. here:
// test vertex location


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
  deallog << "Left vertex=" << (int)tria.begin_active()->face(0)->boundary_indicator() << std::endl;
  deallog << "Right vertex=" << (int)tria.begin_active()->face(1)->boundary_indicator() << std::endl;

  tria.refine_global (2);

  for (typename Triangulation<1,spacedim>::active_cell_iterator
	 cell = tria.begin_active();
       cell != tria.end(); ++cell)
    {
      deallog << "Cell: " << cell << std::endl;
      deallog << "Left vertex=" << (int)cell->face(0)->boundary_indicator() << std::endl;
      deallog << "Right vertex=" << (int)cell->face(1)->boundary_indicator() << std::endl;
    }
}



int main ()
{
  std::ofstream logfile("vertex_as_face_04/output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  test<1> ();
  test<2> ();

  return 0;
}
