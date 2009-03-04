//----------------------------  hp_hanging_node_constraints_02.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2005, 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  hp_hanging_node_constraints_02.cc  ---------------------------


// compute hanging node constraints (before and after the processing
// we do in ConstraintMatrix::close()) for some grids with and without
// random distribution of FEs

#include "../tests.h"
#include <base/logstream.h>
#include <grid/tria.h>
#include <dofs/dof_handler.h>
#include <lac/constraint_matrix.h>
#include <grid/grid_generator.h>

#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_accessor.h>

#include <fe/fe_q.h>

#include <dofs/dof_tools.h>

#include <hp/fe_values.h>
#include <base/quadrature_lib.h>

#include <base/function.h>
#include <numerics/vectors.h>
#include <numerics/matrices.h>

#include <lac/vector.h>
#include <lac/full_matrix.h>
#include <lac/sparse_matrix.h>
#include <lac/solver_cg.h>
#include <lac/precondition.h>

#include <numerics/data_out.h>
#include <fstream>
#include <iostream>

std::ofstream logfile("hp_hanging_nodes_02/output");


template <int dim>
void run (bool random_p,
	  unsigned int *indx) 
{
  Triangulation<dim>     triangulation;
  hp::FECollection<dim>              fe;
  hp::DoFHandler<dim>        dof_handler(triangulation);
  ConstraintMatrix     hanging_node_constraints;

  FE_Q<dim> fe_1 (indx[0]),
    fe_2 (indx[1]),
    fe_3 (indx[2]),
    fe_4 (indx[3]);

  fe.push_back (fe_1);
  fe.push_back (fe_2);
  fe.push_back (fe_3);
  fe.push_back (fe_4);
  
  GridGenerator::hyper_cube (triangulation, -1, 1);
  triangulation.refine_global (5-dim);
  deallog << "Number of active cells: "
	    << triangulation.n_active_cells()
	    << std::endl;
  deallog << "Total number of cells: "
	    << triangulation.n_cells()
	    << std::endl;

				   // Now to the p-Method. Assign
				   // random active_fe_indices to the
				   // different cells.
  typename hp::DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active (),
						     endc = dof_handler.end ();
  if (random_p)
    {      
      for (; cell != endc; ++cell)
	{      
	  cell->set_active_fe_index ((int)(4.0 * (double) random () / (double) RAND_MAX));
	}
    }
  else
    {      
      unsigned int cell_no = 0;  
      for (; cell != endc; ++cell)
	{
	  if (cell_no >= triangulation.n_active_cells () / 2)
	    cell->set_active_fe_index (1);
	  else
	    cell->set_active_fe_index (0);
	}  
    }

  
  dof_handler.distribute_dofs (fe);
  deallog << "Number of degrees of freedom: "
	    << dof_handler.n_dofs()
	    << std::endl;

				   // Create constraints which stem from
				   // the different polynomial degrees on
				   // the different elements.
  hanging_node_constraints.clear ();
  DoFTools::make_hanging_node_constraints (dof_handler,
					   hanging_node_constraints);

  hanging_node_constraints.print (deallog.get_file_stream ());

  hanging_node_constraints.close ();
  hanging_node_constraints.print (deallog.get_file_stream ());
}


template <int dim>
void run_test (unsigned int *indx)
{
  run<dim> (true, indx);
  run<dim> (false, indx);
}



int main () 
{
  logfile.precision(2);
  
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);  

  unsigned int index[] = 
    {
	  1,2,3,4,5,6,7
    };
  
  
  deallog << "Testing Order 1" << std::endl;
  run_test<2> (&(index[0]));
  run_test<3> (&(index[0]));

  deallog << "Testing Order 2" << std::endl;
  run_test<2> (&(index[1]));
  run_test<3> (&(index[1]));

  deallog << "Testing Order 3" << std::endl;
  run_test<2> (&(index[2]));
  run_test<3> (&(index[2]));
      
  return 0;
}
