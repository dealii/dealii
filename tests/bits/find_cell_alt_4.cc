//----------------------------  find_cell_alt_4.cc  ------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2003, 2004, 2005, 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  find_cell_alt_4.cc  ------------------------


// take a 3d mesh and check that we can find an arbitrary point's cell
// in it.

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


void check (Triangulation<3> &tria)
{
   MappingQ1<3> map;
   
  Point<3> p (0.75,0,0);
  
  std::pair<Triangulation<3>::active_cell_iterator, Point<3> > cell
    = GridTools::find_active_cell_around_point (map, tria, p);

  deallog << cell.first << std::endl;
  for (unsigned int v=0; v<GeometryInfo<3>::vertices_per_cell; ++v)
    deallog << "<" << cell.first->vertex(v) << "> ";
  deallog << "[ " << cell.second << "] ";
  deallog << std::endl;

  // Now transform back and check distance
  Point<3> pp = map.transform_unit_to_real_cell(cell.first, GeometryInfo<3>::project_to_unit_cell(cell.second));
  deallog << pp.distance(p) << std::endl;
  Assert (pp.distance(p) < 1e-13,
          ExcInternalError());
}


int main () 
{
  std::ofstream logfile("find_cell_alt_4/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  try
    {
      Triangulation<3> coarse_grid;
      GridGenerator::hyper_cube (coarse_grid);
      coarse_grid.refine_global (3);
      check (coarse_grid);
    }
  catch (const std::exception &exc)
    {
				       // we shouldn't get here...
      deallog << "Caught an error..." << std::endl;
      deallog << exc.what() << std::endl;
    }
}

  
  
