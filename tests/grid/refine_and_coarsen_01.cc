//----------------------------  refine_and_coarsen_3d.cc  ---------------------------
//    $Id: refine_and_coarsen_for_parents_02.cc 23710 2011-05-17 04:50:10Z bangerth $
//    Version: $Name$ 
//
//    Copyright (C) 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-----------------------  refine_and_coarsen_for_parents_02.cc  ----------------------


// check that refine_and_coarsen_fixed_fraction behaves correctly
// if all the indicators are the same.

#include "../tests.h"

#include <deal.II/base/geometry_info.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>

#include <fstream>



template <int dim>
void check ()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube (tria);
  tria.refine_global (2);

  Vector<float> estimated_error_per_cell;
  estimated_error_per_cell.reinit(tria.n_active_cells());
  for (unsigned int j=0; j<estimated_error_per_cell.size(); ++j)
      estimated_error_per_cell(j)= 1.;

  deallog << "n_active_cells: " << tria.n_active_cells() << std::endl;
  
  
  GridRefinement::refine_and_coarsen_fixed_fraction (tria,
						     estimated_error_per_cell,
						     0.25, 0);
  tria.execute_coarsening_and_refinement ();

  deallog << "n_active_cells: " << tria.n_active_cells() << std::endl;
  
  deallog << "OK for " << dim << "d" << std::endl;
}


int main () 
{
  std::ofstream logfile("refine_and_coarsen_01/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  check<1> ();
  check<2> ();
  check<3> ();
}

  
  
