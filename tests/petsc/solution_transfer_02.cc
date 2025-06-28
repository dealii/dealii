// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// SolutionTransfer writes into vectors, but used to forget to
// compress the resulting vectors appropriately. This didn't matter
// most of the time because the (sequential) SolutionTransfer was
// apparently used only with our own vector types, for which
// compress() is a no-op. But it failed with PETSc vectors.
//
// This variation of the function deals with the
// SolutionTransfer::interpolate() variant that takes multiple vectors
// to transfer.


#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/petsc_vector.h>

#include <deal.II/numerics/solution_transfer.h>

#include <iostream>

#include "../tests.h"


template <int dim>
void
test()
{
  Triangulation<dim> triangulation;

  GridGenerator::hyper_cube(triangulation, 0, 1);
  triangulation.refine_global(2);

  FE_Q<dim>       fe(1);
  DoFHandler<dim> dof_handler(triangulation);
  dof_handler.distribute_dofs(fe);


  PETScWrappers::MPI::Vector solution;
  solution.reinit(dof_handler.locally_owned_dofs(), MPI_COMM_SELF);

  // Do a fake transfer
  SolutionTransfer<dim, PETScWrappers::MPI::Vector> solution_trans(dof_handler);

  PETScWrappers::MPI::Vector previous_solution;
  previous_solution = solution;
  triangulation.prepare_coarsening_and_refinement();
  solution_trans.prepare_for_coarsening_and_refinement(previous_solution);

  triangulation.execute_coarsening_and_refinement();

  solution_trans.interpolate(solution);

  deallog << "OK" << std::endl;
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());

  mpi_initlog();
  deallog << std::setprecision(2);

  {
    deallog.push("2d");
    test<2>();
    deallog.pop();
  }

  {
    deallog.push("3d");
    test<3>();
    deallog.pop();
  }
}
