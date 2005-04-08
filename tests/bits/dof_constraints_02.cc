//----------------------------  dof_constraints_02.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2004, 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  dof_constraints_02.cc  ---------------------------


// check DoFConstraints::distribute_local_to_global for vectors

#include "../tests.h"
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

                                   // refine the mesh in a random way so as to
                                   // generate as many hanging node
                                   // constraints as possible
  triangulation.refine_global (4-dim);
  for (unsigned int i=0; i<11-2*dim; ++i)
    {
      typename Triangulation<dim>::active_cell_iterator
        cell = triangulation.begin_active();
      for (unsigned int index=0; cell != triangulation.end(); ++cell, ++index)
        if (index % (3*dim) == 0)
          cell->set_refine_flag();
      triangulation.execute_coarsening_and_refinement ();
    }
  deallog << "Number of cells: " << triangulation.n_active_cells()
          << std::endl;
  
                                   // set up a DoFHandler and compute hanging
                                   // node constraints
  FE_Q<dim> fe(1);
  DoFHandler<dim> dof_handler (triangulation);
  dof_handler.distribute_dofs (fe);
  deallog << "Number of dofs: " << dof_handler.n_dofs() << std::endl;

  ConstraintMatrix constraints;
  DoFTools::make_hanging_node_constraints (dof_handler, constraints);
  constraints.close ();
  deallog << "Number of constraints: " << constraints.n_constraints() << std::endl;

                                   // then set up two vectors
  Vector<double> A(dof_handler.n_dofs()), B(dof_handler.n_dofs());

                                   // then fill the two vectors by setting up
                                   // bogus matrix entries and (1) writing
                                   // them into the vector and condensing away
                                   // hanging node constraints later on, or
                                   // (2) distributing them right away
  std::vector<unsigned int> local_dofs (fe.dofs_per_cell);
  Vector<double> local_vector (fe.dofs_per_cell);
  for (typename DoFHandler<dim>::active_cell_iterator
         cell = dof_handler.begin_active();
       cell != dof_handler.end(); ++cell)
    {
      cell->get_dof_indices (local_dofs);
      local_vector = 0;
      for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
        local_vector(i) = (i+1.)*(local_dofs[i]+1.);

                                       // copy local to global by ourselves
      for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
        A(local_dofs[i]) += local_vector(i);

                                       // or let other functions do that
      constraints.distribute_local_to_global (local_vector, local_dofs, B);
    }

                                   // now condense away constraints from A
  constraints.condense (A);

                                   // we haven't yet set the entries
                                   // for constrained nodes. we can do so at
                                   // will, since these values don't matter
                                   // anyway
  for (unsigned int i=0; i<dof_handler.n_dofs(); ++i)
    if (constraints.is_constrained(i))
      B(i) = A(i);
  
                                   // now comes the check: we subtract B from
                                   // A, and make sure that the result is zero
  A -= B;
  deallog << "|A|=" << A.l2_norm() << std::endl;
  deallog << "|B|=" << B.l2_norm() << std::endl;
  Assert (A.l2_norm() < 1e-12*B.l2_norm(),
          ExcInternalError());
}



int main ()
{
  std::ofstream logfile("dof_constraints_02.output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

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
