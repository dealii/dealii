//----------------------------  constraint_graph.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2005, 2006, 2007, 2008 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  constraint_graph.cc  ---------------------------


// take hp/random and generate the graph of constraints


#include "../tests.h"
#include <base/logstream.h>
#include <grid/tria.h>
#include <grid/grid_generator.h>
#include <grid/tria_accessor.h>
#include <grid/grid_out.h>
#include <grid/tria_iterator.h>
#include <hp/dof_handler.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>
#include <lac/constraint_matrix.h>
#include <fe/fe_q.h>

#include <fstream>


template <int dim>
void test ()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global (1);
  tria.begin_active()->set_refine_flag ();
  tria.execute_coarsening_and_refinement ();
  tria.refine_global (1);  

  hp::FECollection<dim> fe_collection;
  fe_collection.push_back(FE_Q<dim> (1));
  fe_collection.push_back(FE_Q<dim> (2));
  fe_collection.push_back(FE_Q<dim> (3));
  fe_collection.push_back(FE_Q<dim> (4));

  hp::DoFHandler<dim> dof_handler(tria);

  for (typename hp::DoFHandler<dim>::active_cell_iterator
	 cell = dof_handler.begin_active();
       cell != dof_handler.end(); ++cell)
    cell->set_active_fe_index (rand() % fe_collection.size());
  
  dof_handler.distribute_dofs(fe_collection);

  ConstraintMatrix cm;
  DoFTools::make_hanging_node_constraints (dof_handler, cm);
  cm.write_dot (deallog.get_file_stream());
}


int main ()
{
  std::ofstream logfile("constraint_graph/output");
  logfile << std::setprecision(2);
  deallog << std::setprecision(2);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);  

  test<3> ();
  
  deallog << "OK" << std::endl;
}
