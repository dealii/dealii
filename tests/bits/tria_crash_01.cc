//----------------------------  tria_crash_01.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  tria_crash_01.cc  ---------------------------

// a test that checks for a crash introduced in the triangulation class in the
// last few days when fixing refine_and_coarsen_3d


#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>

#include <fstream>
#include <iomanip>


bool predicate (const Point<3> &p,
		const double    diameter)
{
  return ((p[0]-.2)*(p[0]-.2) + (p[2]-p[1]/4)*(p[2]-p[1]/4) < diameter * diameter);
}


int main ()
{
  std::ofstream logfile("tria_crash_01/output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  const unsigned int dim=3;
  Triangulation<dim> tria;
  GridGenerator::cylinder(tria, 1, .7);

  deallog << "n_cells=" << tria.n_active_cells() << std::endl;

  tria.refine_global(2);

  deallog << "n_cells=" << tria.n_active_cells() << std::endl;
  
				   // build up a map of vertex indices
				   // of boundary vertices to the new
				   // boundary points
  std::map<unsigned int,Point<dim> > new_points;
  
  Triangulation<dim>::active_cell_iterator cell=tria.begin_active(),
					   endc=tria.end();

  for (cell=tria.begin_active(); cell!=endc; ++cell)
    if (predicate(cell->center(), cell->diameter()))
      cell->set_refine_flag ();
  tria.execute_coarsening_and_refinement();

  deallog << "n_cells=" << tria.n_active_cells() << std::endl;


  for (cell=tria.begin_active(); cell!=endc; ++cell)
    if (!predicate (cell->center(), cell->diameter()))
      cell->set_coarsen_flag ();

				   // make sure there really are no refinement
				   // flags set
  tria.prepare_coarsening_and_refinement();
  for (cell=tria.begin_active(); cell!=endc; ++cell)
    Assert (!cell->refine_flag_set(), ExcInternalError());

  tria.execute_coarsening_and_refinement();

  deallog << "n_cells=" << tria.n_active_cells() << std::endl;
  
  return 0;
}
