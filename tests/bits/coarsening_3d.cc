//----------------------------  coarsening_3d.cc  ---------------------------
//    $Id$
//    Version: 
//
//    Copyright (C) 2003, 2004 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  coarsening_3d.cc  ---------------------------


// this test failed with an internal error somewhere in the coarsening
// functions

#include "../tests.h"
#include "mesh_3d.h"

#include <base/logstream.h>
#include <grid/tria.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/grid_reordering.h>
#include <grid/grid_generator.h>

#include <fstream>


int main () 
{
  std::ofstream logfile("coarsening_3d.output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  Triangulation<3> coarse_grid;
  create_L_shape (coarse_grid);

                                   // refine once, then unrefine again
  coarse_grid.refine_global (1);
  for (Triangulation<3>::active_cell_iterator c=coarse_grid.begin_active();
       c != coarse_grid.end(); ++c)
    c->set_coarsen_flag ();
  coarse_grid.execute_coarsening_and_refinement ();

  deallog << "ok." << std::endl;
}

  
  
