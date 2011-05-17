//----------------------------  refine_and_coarsen_3d.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  refine_and_coarsen_3d.cc  ---------------------------


// check that if we take an locally refined mesh, refine it globally once,
// then coarsen it globally again, that we get the same mesh
//
// Triangulation::fix_coarsen_flags used to be too conservative in allowing
// cells to be coarsened (see today's changes.html entry). in contrast to the
// 1d and 2d cases, the 3d case appears to have an additional problem
// somewhere

#include "../tests.h"

#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>

#include <fstream>



template <int dim>
void check ()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube (tria);
  tria.refine_global (2);
  tria.begin_active()->set_refine_flag();
  tria.execute_coarsening_and_refinement ();

				   // store which cells we have here
  std::vector<typename Triangulation<dim>::active_cell_iterator> cells;
  for (typename Triangulation<dim>::active_cell_iterator cell=tria.begin_active();
       cell != tria.end(); ++cell)
    cells.push_back (cell);

  const unsigned int n_cells = tria.n_active_cells();
  deallog << n_cells << std::endl;

				   // refine the mesh globally, then coarsen
				   // it again globally
  tria.refine_global (1);

  for (typename Triangulation<dim>::active_cell_iterator cell=tria.begin_active();
       cell != tria.end(); ++cell)
    cell->set_coarsen_flag ();
  tria.execute_coarsening_and_refinement ();


				   // verify that we get the same cells again
  deallog << n_cells << ' ' << tria.n_active_cells()
	  << std::endl;
  
  Assert (tria.n_active_cells() == n_cells,
	  ExcInternalError());

  unsigned int index = 0;
  for (typename Triangulation<dim>::active_cell_iterator cell=tria.begin_active();
       cell != tria.end(); ++cell, ++index)
    Assert (cells[index] == cell,
	    ExcInternalError());

}


int main () 
{
  std::ofstream logfile("refine_and_coarsen_3d/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  check<3> ();
}

  
  
