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



// FEValues::get_function_* had a problem when using FE_Nothing

#include <deal.II/base/quadrature_lib.h>

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


int
main()
{
  initlog();

  Triangulation<2> triangulation;
  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global(1);


  hp::FECollection<2> fe_collection;

  fe_collection.push_back(FESystem<2>(FE_Q<2>(1), 1, FE_Q<2>(1), 1));
  fe_collection.push_back(FESystem<2>(FE_Nothing<2>(), 1, FE_Nothing<2>(), 1));

  DoFHandler<2> dof_handler(triangulation);

  // Assign FE to cells
  DoFHandler<2>::active_cell_iterator cell;
  DoFHandler<2>::active_cell_iterator endc = dof_handler.end();


  cell = dof_handler.begin_active();
  cell->set_active_fe_index(1);
  cell++;
  cell->set_active_fe_index(1);
  cell++;
  cell->set_active_fe_index(0);
  cell++;
  cell->set_active_fe_index(0);


  dof_handler.distribute_dofs(fe_collection);

  // Init solution
  Vector<double> solution(dof_handler.n_dofs());
  solution = 1.0;


  SolutionTransfer<2, Vector<double>> solultion_trans(dof_handler);
  solultion_trans.prepare_for_coarsening_and_refinement(solution);

  triangulation.execute_coarsening_and_refinement();
  dof_handler.distribute_dofs(fe_collection);

  Vector<double> new_solution(dof_handler.n_dofs());
  solultion_trans.interpolate(new_solution);

  hp::QCollection<2> q;
  q.push_back(QMidpoint<2>());
  q.push_back(QMidpoint<2>());
  hp::FEValues<2> x_fe_values(fe_collection, q, update_hessians);
  for (DoFHandler<2>::active_cell_iterator cell = dof_handler.begin_active();

       cell != dof_handler.end();
       ++cell)
    {
      x_fe_values.reinit(cell);
      std::vector<Vector<double>> derivatives(
        q[0].size(), Vector<double>(cell->get_fe().n_components()));

      x_fe_values.get_present_fe_values().get_function_laplacians(new_solution,
                                                                  derivatives);
    }

  // we are good if we made it to here
  deallog << "OK" << std::endl;
}
