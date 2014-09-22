// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2014 by the deal.II authors
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



// SolutionTransfer accessed the active_fe_index of the cell from which to
// interpolate to children, but this active_fe_index is now on an inactive
// cell and may no longer be valid. we need to use the number previously
// stored

#include "../tests.h"
#include <fstream>
#include <sstream>
#include <iostream>

#include <deal.II/base/logstream.h>

#include <deal.II/grid/tria.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/hp/dof_handler.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/hp/fe_collection.h>

#include <deal.II/numerics/solution_transfer.h>

#include <deal.II/numerics/data_out.h>

using namespace dealii;


template <int dim>
void test ()
{
  Triangulation<dim> triangulation;
  GridGenerator::hyper_cube (triangulation);

  hp::FECollection<dim>     fe_collection;
  fe_collection.push_back(FE_Q<dim>(1));
  fe_collection.push_back(FE_Q<dim>(2));

  hp::DoFHandler<dim> dof_handler(triangulation);
  dof_handler.begin_active()->set_active_fe_index(0);
  dof_handler.distribute_dofs (fe_collection);

  // Init solution
  Vector<double> solution(dof_handler.n_dofs());
  solution = 1.0;


  // set refine flag for the only cell we have, then do the refinement
  SolutionTransfer<dim, Vector<double>, hp::DoFHandler<dim> >
    solution_trans(dof_handler);
  dof_handler.begin_active()->set_refine_flag ();
  solution_trans.prepare_for_coarsening_and_refinement(solution);
  triangulation.execute_coarsening_and_refinement ();

  // now set the active_fe_index flags on the new set of fine level cells
  for (unsigned int c=0; c<dof_handler.begin(0)->n_children(); ++c)
    dof_handler.begin(0)->child(c)->set_active_fe_index(1);

  // distribute dofs and transfer solution there
  dof_handler.distribute_dofs (fe_collection);
  
  Vector<double> new_solution(dof_handler.n_dofs());
  solution_trans.interpolate(solution, new_solution);

  // we should now have only 1s in the new_solution vector
  for (unsigned int i=0; i<new_solution.size(); ++i)
    Assert (new_solution[i] == 1, ExcInternalError());
  
  // we are good if we made it to here
  deallog << "OK" << std::endl;
}


int main()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test<1> ();
  test<2> ();
  test<3> ();
}
