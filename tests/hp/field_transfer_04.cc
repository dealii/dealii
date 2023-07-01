// ---------------------------------------------------------------------
//
// Copyright (C) 2021 - 2022 by the deal.II authors
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

// Check that we can use FieldTransfer when multiple cells are changed from
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

// FieldTransfer with ghosted vector

int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  MPILogInitAll mpi_log;

  parallel::distributed::Triangulation<2> triangulation(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global(2);

  hp::FECollection<2> fe_collection;
  fe_collection.push_back(FE_Q<2>(1));
  fe_collection.push_back(FE_Nothing<2>());

  DoFHandler<2> dof_handler(triangulation);

  // Assign FE_Nothing to half of the domain
  for (auto cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          if (cell->center()[1] < 0.5)
            {
              cell->set_active_fe_index(0);
            }
          else
            {
              cell->set_active_fe_index(1);
            }
        }
    }

  dof_handler.distribute_dofs(fe_collection);

  // Initialize solution
  auto locally_relevant_dofs =
    DoFTools::extract_locally_relevant_dofs(dof_handler);
  LinearAlgebra::distributed::Vector<double> solution(
    dof_handler.locally_owned_dofs(), locally_relevant_dofs, MPI_COMM_WORLD);
  const double old_value = 1.;
  for (unsigned int i = 0; i < solution.local_size(); ++i)
    solution.local_element(i) = old_value;

  {
    std::stringstream ss;
    solution.print(ss);
    deallog << ss.str() << std::endl;
  }

  parallel::distributed::experimental::
    FieldTransfer<2, LinearAlgebra::distributed::Vector<double>>
      field_transfer(dof_handler);
  // Assign FE_Q to all the cells
  for (auto cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          cell->set_future_fe_index(0);
        }
    }

  triangulation.prepare_coarsening_and_refinement();

  solution.update_ghost_values();
  field_transfer.prepare_for_coarsening_and_refinement(solution, 1);

  triangulation.execute_coarsening_and_refinement();

  dof_handler.distribute_dofs(fe_collection);
  const unsigned int dofs_per_cell =
    dof_handler.get_fe_collection().max_dofs_per_cell();

  AffineConstraints<double> affine_constraints;
  const double              new_value = 2.;

  locally_relevant_dofs.clear();
  locally_relevant_dofs = DoFTools::extract_locally_relevant_dofs(dof_handler);
  LinearAlgebra::distributed::Vector<double> new_solution(
    dof_handler.locally_owned_dofs(), locally_relevant_dofs, MPI_COMM_WORLD);
  field_transfer.interpolate(new_value, affine_constraints, new_solution);

  {
    std::stringstream ss;
    new_solution.print(ss);
    deallog << ss.str() << std::endl;
  }
}
