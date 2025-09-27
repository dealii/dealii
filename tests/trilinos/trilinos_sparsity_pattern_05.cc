// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// test print functions of Trilinos sparsity pattern

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/trilinos_sparsity_pattern.h>

#include "../tests.h"

void
test()
{
  const int dim = 2;
  // setup system
  dealii::parallel::distributed::Triangulation<dim> triangulation(
    MPI_COMM_WORLD);

  GridGenerator::hyper_cube(triangulation);

  triangulation.refine_global(2);

  const FE_Q<dim> fe_system(1);
  DoFHandler<dim> dh(triangulation);
  dh.distribute_dofs(fe_system);

  const IndexSet relevant_partitioning =
    DoFTools::extract_locally_relevant_dofs(dh);

  // generate empty constraints
  AffineConstraints<double> constraints;

  // generate sparsity pattern
  TrilinosWrappers::SparsityPattern sp(dh.locally_owned_dofs(),
                                       dh.locally_owned_dofs(),
                                       relevant_partitioning,
                                       MPI_COMM_WORLD);

  DoFTools::make_sparsity_pattern(dh,
                                  sp,
                                  constraints,
                                  true,
                                  Utilities::MPI::this_mpi_process(
                                    MPI_COMM_WORLD));
  sp.compress();

  // output
  sp.print_gnuplot(deallog.get_file_stream());
  sp.print(deallog.get_file_stream());
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());
  MPILogInitAll log;

  test();
}
