// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// like _07 but do it the other way around: start with a fine mesh, then
// coarsen

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/hp/fe_collection.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/solution_transfer.h>

#include <iostream>
#include <sstream>

#include "../tests.h"



template <int dim>
void
test()
{
  Triangulation<dim> triangulation;
  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global(1);

  hp::FECollection<dim> fe_collection;
  fe_collection.push_back(FE_Q<dim>(1));
  fe_collection.push_back(FE_Q<dim>(2));

  DoFHandler<dim> dof_handler(triangulation);
  for (unsigned int c = 0; c < dof_handler.begin(0)->n_children(); ++c)
    dof_handler.begin(0)->child(c)->set_active_fe_index(1);
  dof_handler.distribute_dofs(fe_collection);

  // Init solution
  Vector<double> solution(dof_handler.n_dofs());
  solution = 1.0;

  // coarsen everything away
  SolutionTransfer<dim, Vector<double>> solution_trans(dof_handler);
  for (unsigned int c = 0; c < dof_handler.begin(0)->n_children(); ++c)
    dof_handler.begin(0)->child(c)->set_coarsen_flag();

  // now set the active_fe_index flags on the only remaining cell
  dof_handler.begin_active()->set_future_fe_index(0);

  solution_trans.prepare_for_coarsening_and_refinement(solution);
  triangulation.execute_coarsening_and_refinement();

  // distribute dofs and transfer solution there
  dof_handler.distribute_dofs(fe_collection);

  Vector<double> new_solution(dof_handler.n_dofs());
  solution_trans.interpolate(new_solution);

  // we should now have only 1s in the new_solution vector
  for (unsigned int i = 0; i < new_solution.size(); ++i)
    AssertThrow(new_solution[i] == 1, ExcInternalError());

  // we are good if we made it to here
  deallog << "OK" << std::endl;
}


int
main()
{
  initlog();

  test<1>();
  test<2>();
  test<3>();
}
