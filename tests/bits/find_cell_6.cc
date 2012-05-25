//----------------------------  find_cell_5.cc  ---------------------------
//    $Id: find_cell_5.cc 25229 2012-03-09 18:34:59Z pauletti $
//    Version: $Name$ 
//
//    Copyright (C) 2003, 2004, 2005, 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  find_cell_5.cc  ---------------------------

// find_active_cell_around_point() uses the closest vertex
// to find the corresponding cell. This fail(ed) if it is
// a hanging node for that cell.
// This test fortunately does not show this error (yet).

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_boundary_lib.h>

#include <fstream>

#include <deal.II/fe/mapping_q1.h>




void check (Triangulation<2> &tria)
{
  Point<2> p (0.25,0.6);
  
  Triangulation<2>::active_cell_iterator cell
    = GridTools::find_active_cell_around_point (tria, p);

  deallog << cell << std::endl;
  for (unsigned int v=0; v<GeometryInfo<2>::vertices_per_cell; ++v)
    deallog << "<" << cell->vertex(v) << "> ";
  deallog << std::endl;

}


int main () 
{
  std::ofstream logfile("find_cell_6/output");
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

  
  
