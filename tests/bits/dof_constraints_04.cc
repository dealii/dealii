//----------------------------  dof_constraints_04.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2004 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  dof_constraints_04.cc  ---------------------------


// simply check what happens when condensing vectors. This test was
// written when I changed a few things in the algorithm

#include "../tests.h"
#include <lac/sparsity_pattern.h>
#include <lac/sparse_matrix.h>
#include <lac/vector.h>
#include <grid/tria.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/grid_tools.h>
#include <grid/grid_generator.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>
#include <dofs/dof_constraints.h>
#include <fe/fe_q.h>
#include <fstream>
#include <iostream>


template <int dim>
void test ()
{
  deallog << dim << "D" << std::endl;
  
  Triangulation<dim> triangulation;
  GridGenerator::hyper_cube (triangulation);

                                   // refine once, then refine first cell to
                                   // create hanging nodes
  triangulation.refine_global (1);
  triangulation.begin_active()->set_refine_flag ();
  triangulation.execute_coarsening_and_refinement ();
  deallog << "Number of cells: " << triangulation.n_active_cells() << std::endl;
  
                                   // set up a DoFHandler and compute hanging
                                   // node constraints for a Q2 element
  FE_Q<dim> fe(2);
  DoFHandler<dim> dof_handler (triangulation);
  dof_handler.distribute_dofs (fe);
  deallog << "Number of dofs: " << dof_handler.n_dofs() << std::endl;

  ConstraintMatrix constraints;
  DoFTools::make_hanging_node_constraints (dof_handler, constraints);
  constraints.close ();
  deallog << "Number of constraints: " << constraints.n_constraints() << std::endl;

  Vector<double> b(dof_handler.n_dofs());
  for (unsigned int i=0; i<dof_handler.n_dofs(); ++i)
    b(i) = (1.+1.*i*i)/3;

                                   // now condense away constraints
  constraints.condense (b);

                                   // and output what we have
  for (Vector<double>::const_iterator i=b.begin(); i!=b.end(); ++i)
    deallog << *i << std::endl;

                                   // now also make sure that the elements in
                                   // constrained rows are zero
  for (unsigned int i=0; i<b.size(); ++i)
    if (constraints.is_constrained(i))
      Assert (b(i) == 0, ExcInternalError());
}



int main ()
{
  std::ofstream logfile("dof_constraints_04.output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  try
    {
      test<1> ();
      test<2> ();
      test<3> ();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      std::cerr << "Exception on processing: " << std::endl
		<< exc.what() << std::endl
		<< "Aborting!" << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      
      return 1;
    }
  catch (...) 
    {
      std::cerr << std::endl << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      std::cerr << "Unknown exception!" << std::endl
		<< "Aborting!" << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      return 1;
    };
}
