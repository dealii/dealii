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



// testcase by Minh Do-Quang: a case where SolutionTransfer got into trouble
// in a couple of places when using FE_Nothing and FESystem.

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

  // Define compression level for output data
  DataOutBase::VtkFlags vtk_flags;
  vtk_flags.compression_level = DataOutBase::CompressionLevel::best_compression;

  // a follow-up error to the one fixed with _04 was that DataOut also got
  // itself confused
  DataOut<2> data_out2;
  data_out2.attach_dof_handler(dof_handler);
  data_out2.add_data_vector(new_solution, "Solution");
  data_out2.build_patches();
  data_out2.set_flags(vtk_flags);
  data_out2.write_vtu(deallog.get_file_stream());

  // we are good if we made it to here
  deallog << "OK" << std::endl;
}
