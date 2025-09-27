// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test to check if SolutionTransfer works in parallel.
// This tests is based off of solution_transfer_11.cc

#include <deal.II/base/function.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/trilinos_vector.h>

#include <deal.II/numerics/solution_transfer.h>

#include <iostream>
#include <vector>

#include "../tests.h"


template <int dim>
void
transfer(const MPI_Comm mpi_communicator)
{
  const unsigned int this_mpi_process =
    Utilities::MPI::this_mpi_process(mpi_communicator);

  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(1);
  GridTools::partition_triangulation(Utilities::MPI::n_mpi_processes(
                                       mpi_communicator),
                                     tria,
                                     SparsityTools::Partitioner::zoltan);

  hp::FECollection<dim> fe;
  fe.push_back(FE_Q<dim>(1));
  fe.push_back(FE_Q<dim>(2));

  DoFHandler<dim> dof_handler(tria);
  dof_handler.begin(0)->child(0)->set_active_fe_index(1);

  TrilinosWrappers::MPI::Vector solution;

  dof_handler.distribute_dofs(fe);
  IndexSet locally_owned_dofs =
    DoFTools::locally_owned_dofs_per_subdomain(dof_handler)[this_mpi_process];
  IndexSet locally_relevant_dofs =
    DoFTools::locally_relevant_dofs_per_subdomain(
      dof_handler)[this_mpi_process];
  solution.reinit(locally_owned_dofs, mpi_communicator);

  for (unsigned int i = 0; i < solution.size(); ++i)
    if (locally_owned_dofs.is_element(i))
      solution(i) = i;

  SolutionTransfer<dim, TrilinosWrappers::MPI::Vector> soltrans(dof_handler);



  typename Triangulation<dim>::active_cell_iterator cell = tria.begin_active(),
                                                    endc = tria.end();
  ++cell;
  ++cell;
  for (; cell != endc; ++cell)
    cell->set_refine_flag();

  TrilinosWrappers::MPI::Vector old_solution;
  old_solution.reinit(locally_owned_dofs,
                      locally_relevant_dofs,
                      mpi_communicator);
  old_solution = solution;

  tria.prepare_coarsening_and_refinement();
  soltrans.prepare_for_coarsening_and_refinement(old_solution);
  tria.execute_coarsening_and_refinement();
  dof_handler.distribute_dofs(fe);
  locally_owned_dofs =
    DoFTools::locally_owned_dofs_per_subdomain(dof_handler)[this_mpi_process];
  solution.reinit(locally_owned_dofs, mpi_communicator);
  soltrans.interpolate(solution);
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;
  const MPI_Comm                   mpi_communicator = MPI_COMM_WORLD;

  deallog << "   1D solution transfer" << std::endl;
  transfer<1>(mpi_communicator);

  deallog << "   2D solution transfer" << std::endl;
  transfer<2>(mpi_communicator);

  deallog << "   3D solution transfer" << std::endl;
  transfer<3>(mpi_communicator);
}
