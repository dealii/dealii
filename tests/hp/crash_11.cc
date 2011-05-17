//----------------------------  crash_11.cc  ---------------------------
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
//----------------------------  crash_11.cc  ---------------------------


// a test where a degree of freedom was constrained multiple times,
// but with different weights. see the hp paper for more on this

char logname[] = "crash_11/output";


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
#include <deal.II/numerics/vectors.h>
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


				   // create a mesh like this (viewed
				   // from top):
				   //
				   // *---*---*
				   // | 2 | 3 |
				   // *---*---*
				   // | 0 | 1 |
				   // *---*---*
  Triangulation<3>     triangulation;
  std::vector<unsigned int> subdivisions (3, 2);
  subdivisions[2] = 1;
  GridGenerator::subdivided_hyper_rectangle (triangulation, subdivisions,
                                             Point<3>(), Point<3>(2,2,1));

  hp::FECollection<3> fe;
  fe.push_back (FE_Q<3>(1));
  fe.push_back (FE_Q<3>(2));
  fe.push_back (FE_Q<3>(3));
  fe.push_back (FE_Q<3>(4));

  hp::DoFHandler<3>        dof_handler(triangulation);

				   // assign polynomial degrees like this:
                                   //
				   // *---*---*
				   // | 1 | 2 |
				   // *---*---*
				   // | 4 | 3 |
				   // *---*---*
                                   //
  hp::DoFHandler<3>::active_cell_iterator
    cell = dof_handler.begin_active();
  cell->set_active_fe_index (0);
  ++cell;
  cell->set_active_fe_index (1);
  ++cell;
  cell->set_active_fe_index (2);
  ++cell;
  cell->set_active_fe_index (3);

  dof_handler.distribute_dofs (fe);
  
                                   // for illustrative purposes, print
                                   // out the numbers of the dofs that
                                   // belong to the shared edge
                                   // (that's the one that has four
                                   // different fe indices associated
                                   // with it)
  for (hp::DoFHandler<3>::active_line_iterator line = dof_handler.begin_active_line();
       line != dof_handler.end_line(); ++line)
    if (line->n_active_fe_indices() == 4)
      {
	deallog << "Shared line: " << line << std::endl;
	for (unsigned int i=0; i<4; ++i)
	  {
	    deallog << "DoF indices for fe_index=" << i << ": ";
	    std::vector<unsigned int> line_dofs (fe[i].dofs_per_line + 2*fe[i].dofs_per_vertex);
	    line->get_dof_indices (line_dofs, i);
	    for (unsigned int j=0; j<fe[i].dofs_per_line + 2*fe[i].dofs_per_vertex; ++j)
	      deallog << line_dofs[j] << ' ';
	    deallog << std::endl;
	  }
      }
  ConstraintMatrix constraints;
  DoFTools::make_hanging_node_constraints (dof_handler,
					   constraints);
  constraints.close ();

  constraints.print (deallog.get_file_stream());
}

