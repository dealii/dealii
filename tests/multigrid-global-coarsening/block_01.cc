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


/**
 * Test global-coarsening multigrid for block vectors.
 */

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <deal.II/multigrid/mg_constrained_dofs.h>
#include <deal.II/multigrid/mg_transfer_global_coarsening.h>

#include "../tests.h"

template <int dim, typename Number = double>
void
test(const unsigned int n_refinements, const unsigned int fe_degree)
{
  using VectorType = LinearAlgebra::distributed::Vector<Number>;

  const unsigned int min_level = 0;
  const unsigned int max_level = n_refinements;

  parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);

  GridGenerator::hyper_cube(tria);
  tria.refine_global(n_refinements);

  const auto trias =
    MGTransferGlobalCoarseningTools::create_geometric_coarsening_sequence(tria);

  MGLevelObject<DoFHandler<dim>> dof_handlers(min_level, max_level);
  MGLevelObject<DoFHandler<dim>> dof_handlers_system(min_level, max_level);
  MGLevelObject<std::shared_ptr<const Utilities::MPI::Partitioner>>
    partitioners(min_level, max_level);
  MGLevelObject<std::shared_ptr<const Utilities::MPI::Partitioner>>
    partitioners_system(min_level, max_level);
  MGLevelObject<MGTwoLevelTransfer<dim, VectorType>> transfers(min_level,
                                                               max_level);
  MGLevelObject<MGTwoLevelTransfer<dim, VectorType>> transfers_system(
    min_level, max_level);

  for (auto l = min_level; l <= max_level; ++l)
    {
      dof_handlers[l].reinit(*trias[l]);
      dof_handlers[l].distribute_dofs(FE_Q<dim>(fe_degree));
      partitioners[l] = std::make_shared<const Utilities::MPI::Partitioner>(
        dof_handlers[l].locally_owned_dofs(),
        DoFTools::extract_locally_relevant_dofs(dof_handlers[l]),
        MPI_COMM_WORLD);


      dof_handlers_system[l].reinit(*trias[l]);
      dof_handlers_system[l].distribute_dofs(
        FESystem<dim>(FE_Q<dim>(fe_degree), dim));
      DoFRenumbering::block_wise(dof_handlers_system[l]);
      partitioners_system[l] =
        std::make_shared<const Utilities::MPI::Partitioner>(
          dof_handlers_system[l].locally_owned_dofs(),
          DoFTools::extract_locally_relevant_dofs(dof_handlers_system[l]),
          MPI_COMM_WORLD);
    }

  for (unsigned int l = min_level; l < max_level; ++l)
    {
      transfers[l + 1].reinit(dof_handlers[l + 1], dof_handlers[l]);
      transfers_system[l + 1].reinit(dof_handlers_system[l + 1],
                                     dof_handlers_system[l]);
    }



  MGTransferGlobalCoarsening<dim, VectorType> transfer(
    transfers, [&](const auto l, auto &vec) { vec.reinit(partitioners[l]); });

  MGTransferGlobalCoarsening<dim, VectorType> transfer_system(
    transfers_system,
    [&](const auto l, auto &vec) { vec.reinit(partitioners_system[l]); });

  MGTransferBlockGlobalCoarsening<dim, VectorType> transfer_block(transfer);

  LinearAlgebra::distributed::BlockVector<Number> vec_block(dim);
  for (unsigned int d = 0; d < dim; ++d)
    vec_block.block(d).reinit(partitioners.back());
  MGLevelObject<LinearAlgebra::distributed::BlockVector<Number>>
    vec_levels_block(min_level, max_level);

  LinearAlgebra::distributed::Vector<Number> vec_system(
    partitioners_system.back());
  MGLevelObject<LinearAlgebra::distributed::Vector<Number>> vec_levels_system(
    min_level, max_level);

  for (unsigned int b = 0, c = 0; b < dim; ++b)
    for (unsigned int i = 0; i < dof_handlers.back().n_dofs(); ++i, ++c)
      {
        vec_block.block(b)[i] = c;
        vec_system[c]         = c;
      }

  // test copy_to_mg()
  transfer_block.copy_to_mg(dof_handlers.back(), vec_levels_block, vec_block);
  transfer_system.copy_to_mg(dof_handlers_system.back(),
                             vec_levels_system,
                             vec_system);

  vec_levels_block.back().print(deallog.get_file_stream());
  vec_levels_system.back().print(deallog.get_file_stream());

  // test copy_from_mg()
  transfer_block.copy_from_mg(dof_handlers.back(), vec_block, vec_levels_block);
  transfer_system.copy_from_mg(dof_handlers_system.back(),
                               vec_system,
                               vec_levels_system);

  vec_block.print(deallog.get_file_stream());
  vec_system.print(deallog.get_file_stream());


  // test restrict_and_add()
  for (unsigned int l = max_level; l > min_level; --l)
    {
      transfer_block.restrict_and_add(l,
                                      vec_levels_block[l - 1],
                                      vec_levels_block[l]);
      transfer_system.restrict_and_add(l,
                                       vec_levels_system[l - 1],
                                       vec_levels_system[l]);

      vec_levels_block[l - 1].print(deallog.get_file_stream());
      vec_levels_system[l - 1].print(deallog.get_file_stream());
    }

  // test prolongate()
  for (unsigned int l = min_level; l < max_level; ++l)
    {
      transfer_block.prolongate(l + 1,
                                vec_levels_block[l + 1],
                                vec_levels_block[l]);
      transfer_system.prolongate(l + 1,
                                 vec_levels_system[l + 1],
                                 vec_levels_system[l]);

      vec_levels_block[l + 1].print(deallog.get_file_stream());
      vec_levels_system[l + 1].print(deallog.get_file_stream());
    }
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    all;

  deallog.precision(8);

  test<2>(2, 1);
}
