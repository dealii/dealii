// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Check that we can use FieldTransfer when coarsening cells with non-uniform
// data

#include <deal.II/distributed/field_transfer.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/hp/fe_collection.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"


class TestFunction : public Function<2>
{
public:
  double
  value(const dealii::Point<2> &p, const unsigned int component) const override
  {
    return std::exp(10 - 10 * p[1]);
  }
};

void
output(const DoFHandler<2>                              &dof_handler,
       const LinearAlgebra::distributed::Vector<double> &vector,
       const std::string                                &file_name)
{
  DataOut<2> data_out;
  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(vector, "solution");
  data_out.build_patches();
  data_out.write_vtu_in_parallel(file_name, MPI_COMM_WORLD);
}

int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  MPILogInitAll mpi_log;

  parallel::distributed::Triangulation<2> triangulation(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global(3);

  hp::FECollection<2> fe_collection;
  fe_collection.push_back(FE_Q<2>(1));
  fe_collection.push_back(FE_Nothing<2>());

  DoFHandler<2> dof_handler(triangulation);

  // Assign FE_Nothing to half of the domain
  for (auto cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          if (cell->center()[0] < 0.5)
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

  // Initialize with a non-constant function
  VectorTools::interpolate(dof_handler, TestFunction(), solution);

  {
    if (false)
      output(dof_handler, solution, "result_0.vtu");

    deallog << "solution before: " << std::endl;
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

  // coarsen some cells
  for (auto cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          if (cell->center()[0] < 0.25 && cell->center()[1] < 0.25)
            {
              cell->set_coarsen_flag();
            }
        }
    }

  triangulation.prepare_coarsening_and_refinement();

  field_transfer.prepare_for_coarsening_and_refinement(solution, 1);

  triangulation.execute_coarsening_and_refinement();

  dof_handler.distribute_dofs(fe_collection);

  locally_relevant_dofs = DoFTools::extract_locally_relevant_dofs(dof_handler);
  LinearAlgebra::distributed::Vector<double> new_solution(
    dof_handler.locally_owned_dofs(), locally_relevant_dofs, MPI_COMM_WORLD);

  AffineConstraints<double> affine_constraints(
    dof_handler.locally_owned_dofs(),
    DoFTools::extract_locally_relevant_dofs(dof_handler));
  DoFTools::make_hanging_node_constraints(dof_handler, affine_constraints);
  affine_constraints.close();

  const double new_value = 11.;

  field_transfer.interpolate(new_value, affine_constraints, new_solution);

  affine_constraints.distribute(new_solution);

  {
    if (false)
      output(dof_handler, new_solution, "result_1.vtu");

    deallog << "solution after coarsening and activation: " << std::endl;
    std::stringstream ss;
    new_solution.print(ss);
    deallog << ss.str() << std::endl;
  }
}
