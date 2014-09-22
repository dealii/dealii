// ---------------------------------------------------------------------
//
// Copyright (C) 2002 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------



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
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test<2>();
  test<3>();
  return 0;
}
