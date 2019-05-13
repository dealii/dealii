// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2018 by the deal.II authors
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

  IndexSet relevant_partitioning(dh.n_dofs());
  DoFTools::extract_locally_relevant_dofs(dh, relevant_partitioning);

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
