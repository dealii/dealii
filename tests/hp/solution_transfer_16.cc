// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Verify that we can run SolutionTransfer when coarsening a cell that has
// FE_Nothing assigned.


#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/hp/fe_collection.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/solution_transfer.h>

#include "../tests.h"


int
main()
{
  initlog();

  Triangulation<2> triangulation(Triangulation<2>::none);
  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global(1);

  hp::FECollection<2> fe_collection;
  fe_collection.push_back(FE_Q<2>(1));
  fe_collection.push_back(FE_Nothing<2>());

  DoFHandler<2> dof_handler(triangulation);

  // Assign FE_Nothing to the first cell
  dof_handler.begin_active()->set_active_fe_index(1);

  dof_handler.distribute_dofs(fe_collection);
  deallog << "Initial number of dofs: " << dof_handler.n_dofs() << std::endl;

  // Initialize solution
  Vector<double> solution(dof_handler.n_dofs());
  solution = 10;

  deallog << "Initial Vector:" << std::endl;
  solution.print(deallog.get_file_stream());

  // Coarsen all cells
  for (const auto &cell : dof_handler.active_cell_iterators())
    cell->set_coarsen_flag();

  // Assign FE_Q(1) to all cells
  for (const auto &cell : dof_handler.active_cell_iterators())
    cell->set_future_fe_index(0);

  triangulation.prepare_coarsening_and_refinement();

  // Interpolate solution
  SolutionTransfer<2, Vector<double>> solution_trans(dof_handler);

  solution_trans.prepare_for_coarsening_and_refinement(solution);

  triangulation.execute_coarsening_and_refinement();

  dof_handler.distribute_dofs(fe_collection);
  deallog << "Final number of dofs: " << dof_handler.n_dofs() << std::endl;

  Vector<double> new_solution(dof_handler.n_dofs());
  solution_trans.interpolate(new_solution);

  deallog << "Vector after solution transfer:" << std::endl;
  new_solution.print(deallog.get_file_stream());
}
