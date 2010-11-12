//----------------------------  vertex_as_face_05.cc  ---------------------------
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
//----------------------------  vertex_as_face_05.cc  ---------------------------

// verify that we can do things like cell->face() in 1d as well. here:
// test DoFHandler accessors for the same thing as the _01 test


#include "../tests.h"
#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <grid/tria_accessor.h>
#include <grid/grid_generator.h>

#include <fe/fe_system.h>
#include <fe/fe_q.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_accessor.h>

#include <fstream>


template <int spacedim>
void test ()
{
  Triangulation<1,spacedim> tria;
  GridGenerator::hyper_cube (tria);

  FESystem<1,spacedim> fe(FE_Q<1,spacedim>(2),2,
			  FE_Q<1,spacedim>(1),1);
  DoFHandler<1,spacedim> dof_handler (tria);
  dof_handler.distribute_dofs (fe);

  deallog << "Coarse mesh:" << std::endl;
  deallog << "Left vertex=" << dof_handler.begin_active()->face(0)->vertex(0) << std::endl;
  deallog << "Right vertex=" << dof_handler.begin_active()->face(1)->vertex(0) << std::endl;

  tria.refine_global (2);
  dof_handler.distribute_dofs (fe);

  for (typename Triangulation<1,spacedim>::active_cell_iterator
	 cell = dof_handler.begin_active();
       cell != dof_handler.end(); ++cell)
    {
      deallog << "Cell: " << cell << std::endl;
      deallog << "Left vertex=" << cell->face(0)->vertex(0) << std::endl;
      deallog << "Right vertex=" << cell->face(1)->vertex(0) << std::endl;
    }
}



int main ()
{
  std::ofstream logfile("vertex_as_face_05/output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  test<1> ();
  test<2> ();

  return 0;
}
