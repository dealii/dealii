// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// test AffineConstraints<double>::make_consistent_in_parallel for case
// where constraints need to be combined between different ranks

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_generator.h>

#include <deal.II/hp/fe_collection.h>

#include <deal.II/lac/affine_constraints.h>

#include "../tests.h"

using namespace dealii;

void
test()
{
  static const int dim      = 3;
  static const int spacedim = 3;

  parallel::distributed::Triangulation<dim, spacedim> triangulation(
    MPI_COMM_WORLD);
  DoFHandler<dim, spacedim> dof_handler(triangulation);

  hp::FECollection<dim, spacedim> fe_collection;

  IndexSet locally_active_dofs;
  IndexSet locally_relevant_dofs;

  std::vector<unsigned int> repetitions({4, 2, 1});
  Point<dim>                bottom_left(-2, -1, 0);
  Point<dim>                top_right(2, 1, 1);

  GridGenerator::subdivided_hyper_rectangle(triangulation,
                                            repetitions,
                                            bottom_left,
                                            top_right);

  // hp-refine center part
  for (const auto &cell : dof_handler.active_cell_iterators() |
                            IteratorFilters::LocallyOwnedCell())
    {
      cell->set_active_fe_index(1);

      const auto &center = cell->center();
      if (std::abs(center[0]) < 1.)
        {
          if (center[1] > 0.)
            cell->set_active_fe_index(2);
          else
            cell->set_refine_flag();
        }
    }

  triangulation.execute_coarsening_and_refinement();

  // set dofs
  for (unsigned int degree = 1; degree <= 10; ++degree)
    fe_collection.push_back(FE_Q<dim>(degree));
  dof_handler.distribute_dofs(fe_collection);

  const IndexSet &locally_owned_dofs = dof_handler.locally_owned_dofs();
  locally_active_dofs   = DoFTools::extract_locally_active_dofs(dof_handler);
  locally_relevant_dofs = DoFTools::extract_locally_relevant_dofs(dof_handler);
  AffineConstraints<double> constraints(locally_owned_dofs,
                                        locally_relevant_dofs);

  const auto show_constraints_63_64 = [&]() {
    deallog << "What process "
            << Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)
            << " believes about DoF 63:" << std::endl;
    for (const auto &c : constraints.get_lines())
      if (c.index == 63)
        for (const auto &entry : c.entries)
          deallog << "    constrained against " << entry.first
                  << " with weight " << entry.second << std::endl;
    deallog << "What process "
            << Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)
            << " believes about DoF 64:" << std::endl;
    for (const auto &c : constraints.get_lines())
      if (c.index == 64)
        for (const auto &entry : c.entries)
          deallog << "    constrained against " << entry.first
                  << " with weight " << entry.second << std::endl;
  };


  deallog << "------------- make_hanging_node_constraints():" << std::endl;
  DoFTools::make_hanging_node_constraints(dof_handler, constraints);
  show_constraints_63_64();

  deallog << "------------- make_consistent_in_parallel():" << std::endl;
  constraints.make_consistent_in_parallel(locally_owned_dofs,
                                          locally_relevant_dofs,
                                          dof_handler.get_communicator());
  show_constraints_63_64();
}



int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);
  MPILogInitAll                    all;

  test();
}
