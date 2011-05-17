//  ----------------------- find_cell_alt_6.cc  --------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2006 by the deal.II authors and Ralf B. Schulz
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//  ----------------------- find_cell_alt_6.cc  --------------------------


// On a 2D mesh of the following structure look for the cells surrounding
// each vertex, using the find_active_cell_around_point with Mapping:
//
// x-----x-----x
// |     |     |
// |     |     |
// |     |     |
// x--x--x-----x
// |  |  |     |
// x--x--x     x
// |  |  |     |
// x--x--x-----x

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_boundary_lib.h>

#include <deal.II/fe/mapping_q1.h>

#include <fstream>




void check (Triangulation<2> &tria)
{
   const std::vector<Point<2> > &v = tria.get_vertices();
   MappingQ1<2> map;
   
   for(unsigned i=0; i<tria.n_vertices(); i++)
      {
         std::pair<Triangulation<2>::active_cell_iterator, Point<2> >
            cell = GridTools::find_active_cell_around_point(map, tria, v[i]);

         deallog << "Vertex <" << v[i] << "> found in cell ";
         for (unsigned int v=0; v<GeometryInfo<2>::vertices_per_cell; ++v)
            deallog << "<" << cell.first->vertex(v) << "> ";
         deallog << " [local: " << cell.second << "]" << std::endl;
      }
}


int main () 
{
  std::ofstream logfile("find_cell_alt_6/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  try
    {
      Triangulation<2> coarse_grid;
      GridGenerator::hyper_cube (coarse_grid);
      coarse_grid.refine_global (1);
      coarse_grid.begin_active()->set_refine_flag();
      coarse_grid.execute_coarsening_and_refinement();
      check (coarse_grid);
    }
  catch (const std::exception &exc)
    {
				       // we shouldn't get here...
      deallog << "Caught an error..." << std::endl;
      deallog << exc.what() << std::endl;
    }
}

  
  
