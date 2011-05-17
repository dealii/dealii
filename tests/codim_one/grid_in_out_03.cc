
//----------------------------  template.cc  ---------------------------
//    $Id: grid_in_out_02.cc 22595 2010-11-03 11:30:11Z bangerth $
//    Version: $Name$
//
//    Copyright (C) 2005, 2008, 2010, 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  template.cc  ---------------------------


// like grid_in_out but test persistent triangulation 

#include "../tests.h"
#include <fstream>
#include <deal.II/base/logstream.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/persistent_tria.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>

#include <string>

std::ofstream logfile("grid_in_out_03/output");

template <int dim, int spacedim>
void test(std::string filename) {
  logfile << " Tria<" << dim << "," << spacedim << ">: "
	  << filename << std::endl;
  Triangulation<dim, spacedim> tria;
  GridIn<dim, spacedim> gi;
  gi.attach_triangulation (tria);
  std::ifstream in (filename.c_str());
  gi.read_ucd (in);
  
  PersistentTriangulation<dim, spacedim> ptria(tria);
  typename Triangulation<dim,spacedim>::active_cell_iterator cell;

  ptria.restore();
  for(unsigned int i=0; i<2; ++i) 
    {
      for(cell=ptria.begin_active(); cell != ptria.end(); ++cell)
	if(cell->center()[0]<.5)
	  cell->set_refine_flag();
      ptria.execute_coarsening_and_refinement();
    }
  
  ptria.write_flags(logfile);
}

int main ()
{
  deallog.attach(logfile);
  deallog.depth_console(0);

  test<2,3>("grids/square.inp");
  test<2,3>("grids/sphere_1.inp");

  return 0;
}

