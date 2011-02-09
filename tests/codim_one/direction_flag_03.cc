//----------------------------  direction_flag_03.cc  ---------------------------
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
//----------------------------  direction_flag_03.cc  ---------------------------

// there are no direction flags if dim==spacedim. make sure we always
// get back true in that case

#include "../tests.h"

#include <grid/tria.h>
#include <grid/grid_generator.h>
#include <grid/tria_boundary_lib.h>
#include <grid/grid_out.h>
#include <grid/grid_tools.h>

using namespace std;


template <int dim>
void test ()
{
  Triangulation<dim> volume_mesh;
  GridGenerator::hyper_cube(volume_mesh);
  volume_mesh.refine_global(1);

  for (typename Triangulation<dim>::active_cell_iterator
	 cell = volume_mesh.begin_active();
       cell != volume_mesh.end(); ++cell)
    Assert (cell->direction_flag() == true,
	    ExcInternalError());

  deallog << "OK" << std::endl;
}



int main ()
{
  ofstream logfile("direction_flag_03/output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  test<1> ();
  test<2> ();
  test<3> ();
}
