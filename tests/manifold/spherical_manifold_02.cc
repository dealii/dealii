//----------------------------  spherical_manifold_01.cc  ---------------------------
//    Copyright (C) 2011 - 2015 by the mathLab team.
//
//    This file is subject to LGPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  spherical_manifold_02.cc  ---------------------------


// Test that the flat manifold does what it should on a sphere. 

#include "../tests.h"

#include <fstream>
#include <deal.II/base/logstream.h>


// all include files you need here
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/grid_out.h>

// Helper function
template <int dim, int spacedim>
void test(unsigned int ref=1)
{
  SphericalManifold<dim,spacedim> manifold;
  
  Triangulation<dim,spacedim> tria;
  GridGenerator::hyper_ball (tria);
  
  typename Triangulation<dim,spacedim>::active_cell_iterator cell;
  
  for(cell = tria.begin_active(); cell != tria.end(); ++cell) 
    cell->set_all_manifold_ids(1);
  
  for(cell = tria.begin_active(); cell != tria.end(); ++cell) {
    if(cell->center().distance(Point<spacedim>()) < 1e-10)
      cell->set_all_manifold_ids(0);
  }
  
  tria.set_manifold(1, manifold);
  tria.refine_global(2);
  
  GridOut gridout;
  gridout.write_msh(tria, deallog.get_file_stream());
}

int main ()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);
  
  test<2,2>();
  test<3,3>();
  
  return 0;
}

