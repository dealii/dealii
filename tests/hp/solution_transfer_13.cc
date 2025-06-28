// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Verify that we can run SolutionTransfer with FE_Nothing. This led
// to an exception because at one point we try to multiply a vector by
// an empty matrix (because FE_Nothing has no degrees of freedom)
//
// Test case by Claire Bruna-Rosso based on an earlier one written by
// K. Bzowski


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

#include <iostream>

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
  // fe_collection.push_back(FE_Q<2>(2));
  fe_collection.push_back(FE_Nothing<2>());

  DoFHandler<2> dof_handler(triangulation);

  // Assign FE
  DoFHandler<2>::active_cell_iterator cell = dof_handler.begin_active();
  DoFHandler<2>::active_cell_iterator endc = dof_handler.end();

  /*
   * -----------
   * |  0 |  0 |
   * -----------
   * |  1 |  1 |          0 - FEQ, 1 - FE_Nothing
   * -----------
   */

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


  // Vector to visualize the FE of each cell
  Vector<double> FE_Type(triangulation.n_active_cells());
  unsigned int   cnt_cells(0);
  cell = dof_handler.begin_active(), endc = dof_handler.end();
  for (; cell != endc; ++cell)
    {
      unsigned int fe_index = cell->active_fe_index();
      FE_Type[cnt_cells]    = fe_index;
      ++cnt_cells;
    }

  // Save output
  DataOut<2> data_out;
  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(solution, "Solution");
  data_out.add_data_vector(FE_Type, "FE_Type");
  data_out.build_patches();
  data_out.write_gnuplot(deallog.get_file_stream());


  /* Set refine flags:
   * -----------
   * |    |  R |
   * -----------
   * |  R |    |
   * -----------
   */


  cell = dof_handler.begin_active();
  cell->set_refine_flag();
  cell++;
  cell++;
  cell++;
  cell->set_refine_flag();

  triangulation.prepare_coarsening_and_refinement();

  // Interpolate solution
  SolutionTransfer<2, Vector<double>> solution_trans(dof_handler);

  solution_trans.prepare_for_coarsening_and_refinement(solution);

  triangulation.execute_coarsening_and_refinement();

  dof_handler.distribute_dofs(fe_collection);

  Vector<double> new_solution(dof_handler.n_dofs());
  solution_trans.interpolate(new_solution);

  FE_Type.reinit(triangulation.n_active_cells());
  cnt_cells = 0;
  cell = dof_handler.begin_active(), endc = dof_handler.end();
  for (; cell != endc; ++cell)
    {
      unsigned int fe_index = cell->active_fe_index();
      FE_Type[cnt_cells]    = fe_index;
      ++cnt_cells;
    }

  // Save new solution
  DataOut<2> data_out2;
  data_out2.attach_dof_handler(dof_handler);
  data_out2.add_data_vector(new_solution, "Solution");
  data_out2.add_data_vector(FE_Type, "FE_type");
  data_out2.build_patches();
  data_out2.write_gnuplot(deallog.get_file_stream());

  // Solution reinitialization

  dof_handler.distribute_dofs(fe_collection);
  solution.reinit(dof_handler.n_dofs());
  solution = 1.0;

  /* Set coarsen flags:
   * -----------
   * |    |  C |
   * -----------
   * |  C |    |
   * -----------
   */

  endc = dof_handler.end();
  for (cell = dof_handler.begin_active(); cell != endc; ++cell)
    if (cell->level() > 1)
      cell->set_coarsen_flag();

  triangulation.prepare_coarsening_and_refinement();

  // Interpolate solution
  SolutionTransfer<2, Vector<double>> solution_trans2(dof_handler);
  solution_trans2.prepare_for_coarsening_and_refinement(solution);

  triangulation.execute_coarsening_and_refinement();

  dof_handler.distribute_dofs(fe_collection);

  Vector<double> new_solution2(dof_handler.n_dofs());
  solution_trans2.interpolate(new_solution2);

  FE_Type.reinit(triangulation.n_active_cells());
  cnt_cells = 0;
  cell = dof_handler.begin_active(), endc = dof_handler.end();
  for (; cell != endc; ++cell)
    {
      unsigned int fe_index = cell->active_fe_index();
      FE_Type[cnt_cells]    = fe_index;
      ++cnt_cells;
    }

  // Save new solution
  DataOut<2> data_out3;
  data_out3.attach_dof_handler(dof_handler);
  data_out3.add_data_vector(new_solution2, "Solution");
  data_out3.add_data_vector(FE_Type, "FE_type");
  data_out3.build_patches();
  data_out3.write_gnuplot(deallog.get_file_stream());
}
