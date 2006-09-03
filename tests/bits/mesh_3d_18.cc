//----------------------------  mesh_3d_18.cc  ---------------------------
//    $Id: mesh_3d_18.cc 12732 2006-03-28 23:15:45Z wolf $
//    Version: $Name$ 
//
//    Copyright (C) 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  mesh_3d_18.cc  ---------------------------


// adapted from hp/crash_06, which itself is from
// make_hanging_node_constraints for hp elements. used to crash. triggers the
// crash that at the time of writing the test afflicts all
// hp/hp_constraints_*_03 tests

char logname[] = "mesh_3d_18/output";


#include "../tests.h"
#include <base/function.h>
#include <base/logstream.h>

#include <grid/tria.h>
#include <grid/grid_generator.h>
#include <grid/grid_refinement.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>

#include <fstream>
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

  deallog << "Cell=" << cell << std::endl;
  deallog << "Neighbor=" << cell->neighbor(3) << std::endl;
  
  const Triangulation<3>::active_cell_iterator neighbor_child
    = cell->neighbor_child_on_subface (3, 1);

  deallog << "Neighbor_child(3,1)=" << neighbor_child << std::endl;
  for (unsigned int i=0; i<6; ++i)
    {
      deallog << "Neighbor_child neighbors=";
      if (neighbor_child->at_boundary(i))
	deallog << "(at boundary)";
      else
	deallog << neighbor_child->neighbor(i);
      deallog << std::endl;
    }

  deallog << "Neighbor_of_neighbor=" << cell->neighbor_of_neighbor(3) << std::endl;

  deallog << "Face=" << cell->face(3) << std::endl;
  for (unsigned int i=0; i<4; ++i)
    deallog << "Face_child=" << cell->face(3)->child(i) << std::endl;
  
  for (unsigned int i=0; i<4; ++i)
    deallog << "Neighbor_face=" << cell->neighbor(3)->face(i) << std::endl;

  for (unsigned int i=0; i<4; ++i)
    deallog << "Neighbor_child_face=" << neighbor_child->face(i) << std::endl;
  
  Assert (cell->face(3)->child(1) ==
	  neighbor_child->face(cell->neighbor_of_neighbor(3)),
	  ExcInternalError());
}



int main ()
{
  std::ofstream logfile(logname);
  logfile.precision (3);
  
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test_with_wrong_face_orientation ();

  deallog << "OK" << std::endl;
}

