//----------------------------  periodicity_03.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2002, 2003, 2004, 2005, 2010, 2013 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  periodicity_03.cc  ---------------------------


// check periodic boundary conditions for a simple enough case where we know
// the exact set of constraints
//
// this test simply uses two hypercubes, refines one of them twice and matches
// the faces at the far ends. this requires recursing into children

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>

#include <iomanip>
#include <fstream>



template <int dim>
void test () 
{
  deallog << dim << "D" << std::endl;
  
  // create a 2x1 (or 2x1x1) mesh and refine the leftmost cell twice
  Triangulation<dim> triangulation;
  std::vector<unsigned int> repetitions (dim, 1);
  repetitions[0] = 2;
  GridGenerator::subdivided_hyper_rectangle (triangulation,
					     repetitions,
					     Point<dim>(),
					     (dim == 2 ?
					      Point<dim>(2,1) :
					      Point<dim>(2,1,1)));
  triangulation.begin_active()->set_refine_flag ();
  triangulation.execute_coarsening_and_refinement ();
  triangulation.begin_active(1)->set_refine_flag ();
  triangulation.execute_coarsening_and_refinement ();
  
  FE_Q<dim>          fe(1);
  DoFHandler<dim>    dof_handler (triangulation);
  dof_handler.distribute_dofs (fe);

  ConstraintMatrix cm;
  DoFTools::make_periodicity_constraints (dof_handler.begin(0)->face(0),
					  (++dof_handler.begin(0))->face(1),
					  cm);
  cm.print (deallog.get_file_stream());
}

    

int main () 
{
  std::ofstream logfile("periodicity_03/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test<2>();
  test<3>();
  return 0;
}
