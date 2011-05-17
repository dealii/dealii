//----------------------------  find_cell_alt_7.cc  ------------------------
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
//----------------------------  find_cell_alt_7.cc  ------------------------


// Use a circular domain with boundary mappings, and determine cell for points
// close to the boundary

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_boundary_lib.h>

#include <deal.II/fe/mapping_q.h>

#include <fstream>


void check (Triangulation<2> &tria)
{
   MappingQ<2> map(5);

   // Test for a number of points, every ten degrees
   for(unsigned int i=0; i<200; i++)
      {
         Point<2> p(std::sin((double)i/100.*M_PI), std::cos((double)i/100.*M_PI));
         p *= 1.-1e-8;

         std::pair<Triangulation<2>::active_cell_iterator, Point<2> > cell
            = GridTools::find_active_cell_around_point (map, tria, p);

         deallog << cell.first << std::endl;
         for (unsigned int v=0; v<GeometryInfo<2>::vertices_per_cell; ++v)
            deallog << "<" << cell.first->vertex(v) << "> ";
         deallog << "[ " << cell.second << "] ";

         // Now transform back and check distance
         Point<2> pp = map.transform_unit_to_real_cell(cell.first, GeometryInfo<2>::project_to_unit_cell(cell.second));
         deallog << pp.distance(p) << std::endl;
         Assert (pp.distance(p) < 1e-13,
                 ExcInternalError());
      }
  

}


int main () 
{
  std::ofstream logfile("find_cell_alt_7/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  {  
    Triangulation<2> coarse_grid;
    GridGenerator::hyper_ball (coarse_grid);
    static const HyperBallBoundary<2> boundary;
    coarse_grid.set_boundary (0, boundary);
    coarse_grid.refine_global (2);
    check (coarse_grid);
  }
}

  
  
