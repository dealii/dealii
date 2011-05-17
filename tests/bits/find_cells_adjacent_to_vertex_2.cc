//  ---------------- find_cells_adjacent_to_vertex_2.cc  ------------------
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
//  ---------------- find_cells_adjacent_to_vertex_2.cc  ------------------


// Same as the first test, but on a 3D grid of the same structure

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_boundary_lib.h>

#include <fstream>




void check (Triangulation<3> &tria)
{
   for(unsigned i=0; i<tria.n_vertices(); i++)
      {
         std::vector<Triangulation<3>::active_cell_iterator>
            cells = GridTools::find_cells_adjacent_to_vertex(tria, i);

         deallog << "Vertex " << i << " at " << tria.get_vertices()[i] << ": " << cells.size() << " cells" << std::endl;

         for(unsigned c=0; c<cells.size(); c++) {
            for (unsigned int v=0; v<GeometryInfo<3>::vertices_per_cell; ++v)
               deallog << "<" << cells[c]->vertex(v) << "> ";
            deallog << std::endl;
         }
      }
}


int main () 
{
  std::ofstream logfile("find_cells_adjacent_to_vertex_2/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  try
    {
      Triangulation<3> coarse_grid;
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

  
  
