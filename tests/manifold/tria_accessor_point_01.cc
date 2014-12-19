//----------------------------  spherical_manifold_01.cc  ---------------------------
//    Copyright (C) 2011, 2013, 2014 by the mathLab team.
//
//    This file is subject to LGPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  spherical_manifold_01.cc  ---------------------------


// Test spherical manifold on hyper shells.

#include "../tests.h"
#include <fstream>
#include <base/logstream.h>


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
  deallog << "Testing dim " << dim 
	  << ", spacedim " << spacedim << std::endl;

  SphericalManifold<dim,spacedim> manifold;
  
  Triangulation<dim,spacedim> tria;
  GridGenerator::hyper_shell (tria, Point<spacedim>(), .3, .6, 12);

  for(typename Triangulation<dim,spacedim>::active_cell_iterator cell = tria.begin_active(); cell != tria.end(); ++cell) {
    cell->set_all_manifold_ids(1);
  }
  
  tria.set_manifold(1, manifold);
  tria.refine_global(1);

  for(typename Triangulation<dim,spacedim>::active_cell_iterator cell = tria.begin_active(); cell != tria.end(); ++cell) {
    for(unsigned int f=0;f<GeometryInfo<dim>::faces_per_cell; ++f)
      if(cell->face(f)->at_boundary())
	deallog << "Center: " << cell->face(f)->center(true) 
		<< ", Norm: " << cell->face(f)->center(true).norm() << std::endl;
  }
  
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

