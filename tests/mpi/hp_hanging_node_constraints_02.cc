// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// this test failed at some stage:
// DoFTools::make_hanging_node_constraints() has not worked with an
// hp::DoFHandler on ghost cells that neighbor artificial cells.
//
// see: https://github.com/dealii/dealii/issues/9517


#include <deal.II/base/index_set.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_generator.h>

#include <deal.II/hp/fe_collection.h>

#include <deal.II/lac/affine_constraints.h>

#include "../tests.h"


template <int dim>
void
test()
{
  // setup triangulation
  parallel::distributed::Triangulation<dim> triangulation(
    MPI_COMM_WORLD,
    typename Triangulation<dim>::MeshSmoothing(
      Triangulation<dim>::smoothing_on_refinement |
      Triangulation<dim>::smoothing_on_coarsening));
  GridGenerator::hyper_rectangle(
    triangulation,
    (dim == 3 ? Point<dim>(0.0, 0.0, 0.0) : Point<dim>(0.0, 0.0)),
    (dim == 3 ? Point<dim>(1.0, 1.0, 1.0) : Point<dim>(1.0, 1.0)),
    true);
  triangulation.refine_global(3);

  // Refine all the cells with y < 0.5 twice
  for (unsigned int n_ref = 0; n_ref < 2; ++n_ref)
    {
      for (const auto &cell : triangulation.active_cell_iterators())
        if (cell->is_locally_owned() && cell->center()[1] < 0.5)
          cell->set_refine_flag();

      triangulation.prepare_coarsening_and_refinement();
      triangulation.execute_coarsening_and_refinement();
    }

  // setup finite elements
  FESystem<dim> void_fe(FE_Nothing<dim>(), dim);
  FESystem<dim> solid_fe(FE_Q<dim>(1), dim);

  hp::FECollection<dim> fe_collection;
  fe_collection.push_back(void_fe);
  fe_collection.push_back(solid_fe);

  DoFHandler<dim> dof_handler(triangulation);

  // Assign void_fe to all the cells with x < 0.5
  for (const auto &cell : dof_handler.active_cell_iterators() |
                            IteratorFilters::LocallyOwnedCell())
    {
      if (cell->center()[0] < 0.5)
        cell->set_active_fe_index(0);
      else
        cell->set_active_fe_index(1);
    }

  dof_handler.distribute_dofs(fe_collection);

  // make constraints
  const IndexSet locally_relevant_dofs =
    DoFTools::extract_locally_relevant_dofs(dof_handler);

  AffineConstraints<double> hanging_node_constraints;
  hanging_node_constraints.reinit(dof_handler.locally_owned_dofs(),
                                  locally_relevant_dofs);
  DoFTools::make_hanging_node_constraints(dof_handler,
                                          hanging_node_constraints);

  deallog << "OK" << std::endl;
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  MPILogInitAll all;

  deallog.push("2d");
  test<2>();
  deallog.pop();
  MPI_Barrier(MPI_COMM_WORLD);
  deallog.push("3d");
  test<3>();
  deallog.pop();
}
