//----------------------------  mesh_3d_19.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2006, 2008 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  mesh_3d_19.cc  ---------------------------


// an attempt to understand the failure of mesh_3d_18

char logname[] = "mesh_3d_19/output";


#include "../tests.h"
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <fstream>
#include <iomanip>
#include <vector>




void test_with_wrong_face_orientation ()
{
  Triangulation<3>     triangulation;
  GridGenerator::hyper_ball (triangulation);
  triangulation.begin_active()->set_refine_flag ();
  triangulation.execute_coarsening_and_refinement ();

  Triangulation<3>::active_cell_iterator cell = triangulation.begin_active();
  ++cell;
  ++cell;

  deallog << "cell=" << cell << std::endl;
  deallog << "cell->neighbor(3)=" << cell->neighbor(3) << std::endl;
  deallog << "cell->face(3)=" << cell->face(3) << std::endl;

  for (unsigned int i=0; i<6; ++i)
    deallog << "cell->neighbor(3)->face(" << i << ")="
	    << cell->neighbor(3)->face(i) << std::endl;
}



int main ()
{
  std::ofstream logfile(logname);
  deallog << std::setprecision (3);
  
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test_with_wrong_face_orientation ();

  deallog << "OK" << std::endl;
}

