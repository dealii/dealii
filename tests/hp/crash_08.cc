// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2013 by the deal.II authors
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



// the crash_08 testcase discussed in the hp paper. this produces a cyclic
// constraint between degrees of freedom 3->14->17->6->3 with the algorithm
// that is presently in make_hanging_node_constraints

char logname[] = "output";


#include "../tests.h"
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/lac/vector.h>

#include <deal.II/grid/tria.h>
#include <deal.II/hp/dof_handler.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/fe/fe_q.h>

#include <fstream>
#include <vector>



int main ()
{
  std::ofstream logfile(logname);
  logfile.precision (3);

  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);


  // create a mesh like this:
  //
  // *---*---*---*
  // | 6 | 7 | 8 |
  // *---*---*---*
  // | 3 | 4 | 5 |
  // *---*---*---*
  // | 0 | 1 | 2 |
  // *---*---*---*
  Triangulation<2>     triangulation;
  GridGenerator::subdivided_hyper_cube (triangulation, 3);

  hp::FECollection<2> fe;
  fe.push_back (FE_Q<2>(1));
  fe.push_back (FE_Q<2>(2));
  fe.push_back (FE_Q<2>(3));

  hp::DoFHandler<2>        dof_handler(triangulation);

  // subdivide cells 1, 3, 5, 7
  hp::DoFHandler<2>::active_cell_iterator
  cell = dof_handler.begin_active();
  ++cell;
  cell->set_refine_flag ();
  ++cell;
  ++cell;
  cell->set_refine_flag ();
  ++cell;
  ++cell;
  cell->set_refine_flag ();
  ++cell;
  ++cell;
  cell->set_refine_flag ();
  triangulation.execute_coarsening_and_refinement ();

  // now set fe_index as described in the
  // paper
  cell = dof_handler.begin_active();
  cell->set_active_fe_index (0);
  ++cell;
  cell->set_active_fe_index (0);
  ++cell;
  cell->set_active_fe_index (1);
  ++cell;
  cell->set_active_fe_index (0);
  ++cell;
  cell->set_active_fe_index (0);

  // one set of small cells
  ++cell;
  cell->set_active_fe_index (0);
  ++cell;
  cell->set_active_fe_index (0);
  ++cell;
  cell->set_active_fe_index (0);
  ++cell;
  cell->set_active_fe_index (2);

  // one set of small cells
  ++cell;
  cell->set_active_fe_index (0);
  ++cell;
  cell->set_active_fe_index (2);
  ++cell;
  cell->set_active_fe_index (0);
  ++cell;
  cell->set_active_fe_index (0);

  // one set of small cells
  ++cell;
  cell->set_active_fe_index (0);
  ++cell;
  cell->set_active_fe_index (0);
  ++cell;
  cell->set_active_fe_index (2);
  ++cell;
  cell->set_active_fe_index (0);

  // one set of small cells
  ++cell;
  cell->set_active_fe_index (2);
  ++cell;
  cell->set_active_fe_index (0);
  ++cell;
  cell->set_active_fe_index (0);
  ++cell;
  cell->set_active_fe_index (0);

  dof_handler.distribute_dofs (fe);

  ConstraintMatrix constraints;
  DoFTools::make_hanging_node_constraints (dof_handler,
                                           constraints);
  constraints.close ();

  constraints.print (deallog.get_file_stream());
}

