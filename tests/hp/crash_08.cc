//----------------------------  crash_08.cc  ---------------------------
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
//----------------------------  crash_08.cc  ---------------------------


// the crash_08 testcase discussed in the hp paper. this produces a cyclic
// constraint between degrees of freedom 3->14->17->6->3 with the algorithm
// that is presently in make_hanging_node_constraints

char logname[] = "crash_08/output";


#include "../tests.h"
#include <base/function.h>
#include <base/logstream.h>
#include <base/quadrature_lib.h>
#include <lac/vector.h>

#include <grid/tria.h>
#include <hp/dof_handler.h>
#include <lac/constraint_matrix.h>
#include <grid/grid_generator.h>
#include <grid/grid_refinement.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/tria_boundary_lib.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>
#include <numerics/vectors.h>
#include <fe/fe_q.h>

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
  ++cell; ++cell;
  cell->set_refine_flag ();
  ++cell; ++cell;
  cell->set_refine_flag ();
  ++cell; ++cell;
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

