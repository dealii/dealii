// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// This test crashed at some point: We have set and sent active_fe_indices based
// on the refinement flags on the p::d::Triangulation object. However, p4est has
// the last word on deciding which cells will be refined -- and p4est makes use
// of it in the specific scenario provided as a test. A fix has been introduced
// along with this test.


#include <deal.II/base/function.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h>

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_generator.h>

#include <deal.II/hp/fe_collection.h>

#include <deal.II/lac/la_parallel_vector.h>

#include <deal.II/numerics/solution_transfer.h>
#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"


template <int dim>
void
test()
{
  // setup
  parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(tria);
  tria.refine_global(2);

  const unsigned int    max_degree = 6 - dim;
  hp::FECollection<dim> fe_dgq;
  for (unsigned int deg = 1; deg <= max_degree; ++deg)
    fe_dgq.push_back(FE_DGQ<dim>(deg));

  DoFHandler<dim> dgq_dof_handler(tria);

  // randomly assign FEs
  for (const auto &cell : dgq_dof_handler.active_cell_iterators() |
                            IteratorFilters::LocallyOwnedCell())
    cell->set_active_fe_index(Testing::rand() % max_degree);
  dgq_dof_handler.distribute_dofs(fe_dgq);

  // prepare index sets
  IndexSet dgq_locally_owned_dofs = dgq_dof_handler.locally_owned_dofs();
  IndexSet dgq_locally_relevant_dofs =
    DoFTools::extract_locally_relevant_dofs(dgq_dof_handler);
  IndexSet dgq_ghost_dofs = dgq_locally_relevant_dofs;
  dgq_ghost_dofs.subtract_set(dgq_locally_owned_dofs);

  // prepare dof_values
  LinearAlgebra::distributed::Vector<double> dgq_solution;
  dgq_solution.reinit(dgq_locally_owned_dofs, dgq_ghost_dofs, MPI_COMM_WORLD);
  VectorTools::interpolate(dgq_dof_handler,
                           Functions::ZeroFunction<dim>(),
                           dgq_solution);
  dgq_solution.update_ghost_values();

  SolutionTransfer<dim, LinearAlgebra::distributed::Vector<double>>


    dgq_soltrans(dgq_dof_handler);
  dgq_soltrans.prepare_for_coarsening_and_refinement(dgq_solution);

  // refine and transfer
  {
    unsigned int counter = 0;
    for (auto cell = tria.begin_active(); cell != tria.end(); ++cell, ++counter)
      if (cell->is_locally_owned())
        {
          if (counter > ((dim == 2) ? 4 : 8))
            cell->set_coarsen_flag();
          else
            cell->set_refine_flag();
        }
  }

  tria.execute_coarsening_and_refinement();
  dgq_dof_handler.distribute_dofs(fe_dgq);

  // prepare index sets
  dgq_locally_owned_dofs = dgq_dof_handler.locally_owned_dofs();
  dgq_locally_relevant_dofs =
    DoFTools::extract_locally_relevant_dofs(dgq_dof_handler);
  dgq_ghost_dofs = dgq_locally_relevant_dofs;
  dgq_ghost_dofs.subtract_set(dgq_locally_owned_dofs);

  // unpack dof_values
  dgq_solution.reinit(dgq_locally_owned_dofs, dgq_ghost_dofs, MPI_COMM_WORLD);
  dgq_soltrans.interpolate(dgq_solution);

  deallog << "OK" << std::endl;
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;

  deallog.push("2d");
  test<2>();
  deallog.pop();
  deallog.push("3d");
  test<3>();
  deallog.pop();
}
