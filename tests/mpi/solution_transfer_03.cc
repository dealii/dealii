// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


// Test to check if SolutionTransfer works in parallel for block vectors.
// This tests is based off of mpi/solution_transfer_02.cc

#include <deal.II/base/function.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/trilinos_parallel_block_vector.h>

#include <deal.II/numerics/solution_transfer.h>

#include <iostream>
#include <vector>

#include "../tests.h"


template <typename DH>
void
initialize_indexsets(IndexSet &             locally_owned_dofs,
                     IndexSet &             locally_relevant_dofs,
                     std::vector<IndexSet> &locally_owned_partitioning,
                     std::vector<IndexSet> &locally_relevant_partitioning,
                     const DH &             dof_handler,
                     const std::vector<unsigned int> &block_component,
                     const unsigned int               this_mpi_process)
{
  locally_owned_dofs =
    DoFTools::locally_owned_dofs_per_subdomain(dof_handler)[this_mpi_process];
  locally_relevant_dofs = DoFTools::locally_relevant_dofs_per_subdomain(
    dof_handler)[this_mpi_process];

  const unsigned int                         n_blocks = block_component.size();
  const std::vector<types::global_dof_index> dofs_per_block =
    DoFTools::count_dofs_per_fe_block(dof_handler, block_component);

  locally_owned_partitioning.clear();
  locally_relevant_partitioning.clear();
  locally_owned_partitioning.reserve(n_blocks);
  locally_relevant_partitioning.reserve(n_blocks);

  for (unsigned int b = 0; b < n_blocks; ++b)
    {
      const types::global_dof_index idx_begin =
        std::accumulate(dofs_per_block.begin(),
                        std::next(dofs_per_block.begin(), b),
                        0);
      const types::global_dof_index idx_end =
        std::accumulate(dofs_per_block.begin(),
                        std::next(dofs_per_block.begin(), b + 1),
                        0);
      locally_owned_partitioning.push_back(
        locally_owned_dofs.get_view(idx_begin, idx_end));
      locally_relevant_partitioning.push_back(
        locally_relevant_dofs.get_view(idx_begin, idx_end));
    }
}


template <int dim>
void
transfer(const MPI_Comm &mpi_communicator)
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

  const FESystem<dim>             fe(FE_Q<dim>(2), 1, FE_Q<dim>(1), 1);
  const std::vector<unsigned int> block_component({0, 1});

  DoFHandler<dim> dof_handler(tria);

  TrilinosWrappers::MPI::BlockVector solution;

  dof_handler.distribute_dofs(fe);
  DoFRenumbering::component_wise(dof_handler, block_component);

  IndexSet              locally_owned_dofs, locally_relevant_dofs;
  std::vector<IndexSet> locally_owned_partitioning,
    locally_relevant_partitioning;
  initialize_indexsets(locally_owned_dofs,
                       locally_relevant_dofs,
                       locally_owned_partitioning,
                       locally_relevant_partitioning,
                       dof_handler,
                       block_component,
                       this_mpi_process);

  solution.reinit(locally_owned_partitioning, mpi_communicator);

  for (unsigned int i = 0; i < solution.size(); ++i)
    if (locally_owned_dofs.is_element(i))
      solution(i) = i;

  SolutionTransfer<dim, TrilinosWrappers::MPI::BlockVector, DoFHandler<dim>>
    soltrans(dof_handler);

  typename Triangulation<dim>::active_cell_iterator cell = tria.begin_active(),
                                                    endc = tria.end();
  ++cell;
  ++cell;
  for (; cell != endc; ++cell)
    cell->set_refine_flag();

  TrilinosWrappers::MPI::BlockVector old_solution;
  old_solution.reinit(locally_owned_partitioning,
                      locally_relevant_partitioning,
                      mpi_communicator);
  old_solution = solution;

  tria.prepare_coarsening_and_refinement();
  soltrans.prepare_for_pure_refinement();
  tria.execute_coarsening_and_refinement();

  dof_handler.distribute_dofs(fe);
  DoFRenumbering::component_wise(dof_handler, block_component);

  initialize_indexsets(locally_owned_dofs,
                       locally_relevant_dofs,
                       locally_owned_partitioning,
                       locally_relevant_partitioning,
                       dof_handler,
                       block_component,
                       this_mpi_process);

  solution.reinit(locally_owned_partitioning, mpi_communicator);
  soltrans.refine_interpolate(old_solution, solution);
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;
  const MPI_Comm &                 mpi_communicator = MPI_COMM_WORLD;

  deallog << "   1D solution transfer" << std::endl;
  transfer<1>(mpi_communicator);

  deallog << "   2D solution transfer" << std::endl;
  transfer<2>(mpi_communicator);

  deallog << "   3D solution transfer" << std::endl;
  transfer<3>(mpi_communicator);
}
