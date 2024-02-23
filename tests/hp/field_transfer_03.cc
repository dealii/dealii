// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Check that we can use FieldTransfer when the coarsest cells are changed from
// FE_Nothing to FE_Q


#include <deal.II/distributed/field_transfer.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/hp/fe_collection.h>

#include "../tests.h"


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  initlog();

  parallel::distributed::Triangulation<2> triangulation(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global(1);
  triangulation.execute_coarsening_and_refinement();

  for (auto &cell : triangulation.active_cell_iterators())
    {
      if (cell->center()[1] < 0.5)
        cell->set_refine_flag();
    }
  triangulation.execute_coarsening_and_refinement();


  hp::FECollection<2> fe_collection;
  fe_collection.push_back(FE_Q<2>(2));
  fe_collection.push_back(FE_Nothing<2>());

  DoFHandler<2> dof_handler(triangulation);

  // Assign FE_Nothing to the coarsest cells
  for (auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->center()[1] > 0.5)
        cell->set_active_fe_index(1);
    }

  dof_handler.distribute_dofs(fe_collection);
  const unsigned int initial_n_dofs = dof_handler.n_dofs();

  // Initialize solution
  LinearAlgebra::distributed::Vector<double> solution(initial_n_dofs);
  const double                               old_value = 10.;
  for (unsigned int i = 0; i < initial_n_dofs; ++i)
    solution[i] = old_value;

  parallel::distributed::experimental::
    FieldTransfer<2, LinearAlgebra::distributed::Vector<double>>
      field_transfer(dof_handler);
  // Assign FE_Q to all the cells
  for (auto &cell : dof_handler.active_cell_iterators())
    {
      cell->set_future_fe_index(0);
    }

  triangulation.prepare_coarsening_and_refinement();

  field_transfer.prepare_for_coarsening_and_refinement(solution, 1);

  triangulation.execute_coarsening_and_refinement();

  dof_handler.distribute_dofs(fe_collection);
  const unsigned int dofs_per_cell =
    dof_handler.get_fe_collection().max_dofs_per_cell();

  AffineConstraints<double> affine_constraints;
  DoFTools::make_hanging_node_constraints(dof_handler, affine_constraints);
  affine_constraints.close();
  const double new_value = -1.;

  const unsigned int final_n_dofs = dof_handler.n_dofs();
  LinearAlgebra::distributed::Vector<double> new_solution(final_n_dofs);
  field_transfer.interpolate(new_value, affine_constraints, new_solution);

  affine_constraints.distribute(new_solution);

  new_solution.print(deallog.get_file_stream());
}
