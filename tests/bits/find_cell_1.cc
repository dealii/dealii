//----------------------------  find_cell_1.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2003 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  find_cell_1.cc  ---------------------------


// take a 2d mesh and check that we can find an arbitrary point's cell
// in it

#include <base/logstream.h>
#include <grid/tria.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/grid_tools.h>
#include <grid/grid_generator.h>
#include <grid/tria_boundary_lib.h>

#include <fstream>




void check (Triangulation<2> &tria)
{
  Point<2> p (1./3., 1./2.);
  
  Triangulation<2>::active_cell_iterator cell
    = GridTools::find_active_cell_around_point (tria, p);

  deallog << cell << std::endl;
  for (unsigned int v=0; v<GeometryInfo<2>::vertices_per_cell; ++v)
    deallog << "<" << cell->vertex(v) << "> ";
  deallog << std::endl;

  Assert (p.distance (cell->center()) < cell->diameter()/2,
	  ExcInternalError());
}


int main () 
{
  std::ofstream logfile("find_cell_1.output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  {  
    Triangulation<2> coarse_grid;
    GridGenerator::hyper_cube (coarse_grid);
    coarse_grid.refine_global (2);
    check (coarse_grid);
  }
  
  {  
    Triangulation<2> coarse_grid;
    GridGenerator::hyper_ball (coarse_grid);
    static const HyperBallBoundary<2> boundary;
    coarse_grid.set_boundary (0, boundary);
    coarse_grid.refine_global (2);
    check (coarse_grid);
  }
}

  
  
