//----------------------------  manifold_id_01.cc  ---------------------------
//    Copyright (C) 2011, 2013 by the mathLab team.
//
//    This file is subject to LGPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  flat_manifold_01.cc  ---------------------------


// Test that the flat manifold does what it should

#include "../tests.h"
#include <fstream>
#include <base/logstream.h>


// all include files you need here
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/grid_out.h>

// Helper function
template <int dim, int spacedim>
void test(unsigned int ref=1)
{
  deallog << "Testing dim=" << dim
          << ", spacedim="<< spacedim << std::endl;

  Triangulation<dim,spacedim> tria;
  GridGenerator::hyper_cube (tria);
  tria.refine_global(1);

  typename Triangulation<dim,spacedim>::active_cell_iterator 
    cell;
  
  for(cell=tria.begin_active(); cell!=tria.end(); ++cell) {

	  // check that FlatManifold returns the middle of the cell. 
	  deallog << "Cell: " << cell << std::endl;
	  if(cell->get_manifold().get_new_point_on_cell(cell).distance(cell->center()) > 1e-6) {
		  deallog << "Default manifold: " << cell->get_manifold().get_new_point_on_cell(cell) << std::endl;
		  deallog << "Center of cell  : " << cell->center() << std::endl;
	  } else {
		deallog << "OK!" << std::endl;
	  }
   }
}

int main ()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test<1,1>();
  test<1,2>();
  test<2,2>();
  test<2,3>();
  test<3,3>();

  return 0;
}

