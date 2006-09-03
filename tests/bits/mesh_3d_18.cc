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

  Triangulation<3>::active_cell_iterator cell = triangulation.begin_active(),
					 endc = triangulation.end();
  
  for (; cell!=endc; ++cell)
    for (unsigned int face=0; face<GeometryInfo<3>::faces_per_cell; ++face)
      if (cell->face(face)->has_children()) 
	{
					   // so now we've found a face of an
					   // active cell that has
					   // children. that means that there
					   // are hanging nodes
	  for (unsigned int c=0; c<GeometryInfo<3>::subfaces_per_face; ++c)
	    {
	      Triangulation<3>::active_cell_iterator neighbor_child
		= cell->neighbor_child_on_subface (face, c);

					       // some sanity checks
					       // -- particularly
					       // useful if you start
					       // to think about faces
					       // with
					       // face_orientation==false
					       // and whether we
					       // really really have
					       // the right face...
	      Assert (cell->face(face)->child(c) ==
		      neighbor_child->face(cell->neighbor_of_neighbor(face)),
		      ExcInternalError());
	    }
	}
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

