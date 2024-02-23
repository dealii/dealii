// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------
//
// This test was written by Sam Cox.

// Tests the behavior of DoFTools::make_flux_sparsity_pattern (DoFHandler,
//                          SparsityPattern, AffineConstraints<double>, bool,
//                          coupling, flux_coupling, subdomain_id)

#include <deal.II/base/point.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/trilinos_sparsity_pattern.h>

#include <list>
#include <set>
#include <sstream>

#include "../tests.h"


void
test()
{
  unsigned int np   = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
  unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  const int    dim  = 2;
  // Setup system
  dealii::parallel::distributed::Triangulation<dim> triangulation(
    MPI_COMM_WORLD);

  GridGenerator::hyper_rectangle(triangulation,
                                 Point<dim>(0, 0),
                                 Point<dim>(1, 1),
                                 true);

  triangulation.refine_global(1);
  // Extra refinement to generate hanging nodes
  for (typename Triangulation<dim>::active_cell_iterator cell =
         triangulation.begin_active();
       cell != triangulation.end();
       ++cell)
    if (cell->center()[0] > 0.49)
      cell->set_refine_flag();

  triangulation.prepare_coarsening_and_refinement();
  triangulation.execute_coarsening_and_refinement();
  triangulation.repartition();

  const FESystem<dim> fe_system(FE_Q<dim>(2), 1, FE_DGQ<dim>(2), 1);

  DoFHandler<dim> dh(triangulation);

  dh.distribute_dofs(fe_system);

  // Couple the internal DoFs of both finite elements.
  // Only couple the face terms of the discontinuous element.
  Table<2, DoFTools::Coupling> coupling(2, 2);
  Table<2, DoFTools::Coupling> flux_coupling(2, 2);

  for (unsigned int i = 0; i < 2; ++i)
    for (unsigned int j = 0; j < 2; ++j)
      {
        coupling[i][j]      = DoFTools::none;
        flux_coupling[i][j] = DoFTools::none;
      }

  coupling[0][0]      = DoFTools::always;
  coupling[1][1]      = DoFTools::always;
  flux_coupling[1][1] = DoFTools::always;

  const IndexSet relevant_partitioning =
    DoFTools::extract_locally_relevant_dofs(dh);

  // Generate hanging node constraints
  AffineConstraints<double> constraints;
  constraints.clear();
  DoFTools::make_hanging_node_constraints(dh, constraints);
  constraints.close();

  // Generate sparsity pattern
  DynamicSparsityPattern sp(relevant_partitioning);
  DoFTools::make_flux_sparsity_pattern(dh,
                                       sp,
                                       constraints,
                                       false,
                                       coupling,
                                       flux_coupling,
                                       Utilities::MPI::this_mpi_process(
                                         MPI_COMM_WORLD));
  SparsityTools::distribute_sparsity_pattern(sp,
                                             dh.locally_owned_dofs(),
                                             MPI_COMM_WORLD,
                                             relevant_partitioning);

  // Output
  MPI_Barrier(MPI_COMM_WORLD);

  deallog.push(Utilities::int_to_string(myid));

  deallog << "**** proc " << myid << ": \n\n";
  deallog << "Sparsity pattern:" << std::endl;
  sp.print_gnuplot(deallog.get_file_stream());

  MPI_Barrier(MPI_COMM_WORLD);
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());
  const unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  MPILogInitAll      log;

  deallog.push(Utilities::int_to_string(myid));

  test();
  deallog.pop();
}
