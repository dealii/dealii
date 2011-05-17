//----------------------------  max_n_cells_01.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2007 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  max_n_cells_01.cc  ---------------------------

// test the max_n_cells argument to GridRefinement::refine_and_coarsen_fixed_number


#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_out.h>

#include <fstream>
#include <iomanip>



template <int dim>
void test ()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube (tria);
  tria.refine_global(2);

  for (unsigned int cycle=0; cycle<7*(4-dim)*(4-dim); ++cycle)
    {
      deallog << "cycle=" << cycle << ", n_cells=" << tria.n_active_cells() << std::endl;

      Vector<float> criteria (tria.n_active_cells());
      for (unsigned int i=0; i<tria.n_active_cells(); ++i)
	criteria(i) = i;

      GridRefinement::refine_and_coarsen_fixed_number (tria,
						       criteria,
						       1./4, 1./64,
						       10000);
      tria.execute_coarsening_and_refinement();
    }
}



int main ()
{
  std::ofstream logfile("max_n_cells_01/output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  test<1> ();
  test<2> ();
  test<3> ();
}
