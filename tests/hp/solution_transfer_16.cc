// ---------------------------------------------------------------------
//
// Copyright (C) 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------


// Verify that we can run SolutionTransfer when coarsening a cell that has
// FE_Nothing assigned.


#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/hp/dof_handler.h>
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

  hp::DoFHandler<2> dof_handler(triangulation);

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

  triangulation.prepare_coarsening_and_refinement();

  // Interpolate solution
  SolutionTransfer<2, Vector<double>, hp::DoFHandler<2>> solution_trans(
    dof_handler);
  solution_trans.prepare_for_coarsening_and_refinement(solution);

  triangulation.execute_coarsening_and_refinement();

  // Assign FE_Q(1) to all cells
  for (const auto &cell : dof_handler.active_cell_iterators())
    cell->set_active_fe_index(0);

  dof_handler.distribute_dofs(fe_collection);
  deallog << "Final number of dofs: " << dof_handler.n_dofs() << std::endl;

  Vector<double> new_solution(dof_handler.n_dofs());
  new_solution = 1.;
  solution_trans.interpolate(solution, new_solution);

  deallog << "Vector after solution transfer:" << std::endl;
  new_solution.print(deallog.get_file_stream());
}
