// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2023 by the deal.II authors
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
 * Test transfer operator for polynomial coarsening in a hp-context.
 *
 * Example:
 *
 * +-----+-----+      +-----+-----+      +-----+-----+
 * |     |     |      |     |     |      |     |     |
 * |  3  |  4  |      |  1  |  2  |      |  1  |  1  |
 * |     |     |  0   |     |     |  1   |     |     |
 * +-----+-----+  ->  +-----+-----+  ->  +-----+-----+
 * |     |     |      |     |     |      |     |     |
 * |  1  |  2  |      |  1  |  1  |      |  1  |  1  |
 * |     |     |      |     |     |      |     |     |
 * +-----+-----+      +-----+-----+      +-----+-----+
 *
 *                 ... with fe_degree in the cells
 * like mg_transfer_p_02.cc but testing the assembled transfer matrices.
 */

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/mpi.templates.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>

#include <deal.II/hp/fe_collection.h>

#include <deal.II/lac/sparsity_tools.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <deal.II/multigrid/mg_constrained_dofs.h>
#include <deal.II/multigrid/mg_transfer_matrix_free.templates.h>

#include <limits>

#include "mg_transfer_util.h"

using namespace dealii;

template <int dim, typename Number>
void
do_test()
{
  parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);

  // create grid
  GridGenerator::subdivided_hyper_cube(tria, 2);

  // data structures needed on all levels
  DoFHandler<dim>       dof_handler_fine(tria);
  DoFHandler<dim>       dof_handler_coarse(tria);
  hp::FECollection<dim> fe_collection;


  // loop over levels
  for (unsigned int l = 0; l < std::numeric_limits<unsigned int>::max(); ++l)
    {
      deallog.push("level" + std::to_string(l));

      // set FEs on fine level
      if (l == 0)
        {
          unsigned int i = 0;

          for (auto &cell : dof_handler_fine.active_cell_iterators())
            {
              if (cell->is_locally_owned())
                cell->set_active_fe_index(i);
              fe_collection.push_back(FE_Q<dim>(i + 1));
              i++;
            }
        }
      else
        {
          for (auto &cell : dof_handler_fine.active_cell_iterators())
            {
              if (cell->is_locally_owned())
                cell->set_active_fe_index(
                  std::max((cell->active_fe_index() + 1u) / 2 /*bisection*/,
                           1u) -
                  1);
            }
        }

      // set FEs on coarse level
      {
        auto cell_other = dof_handler_fine.begin_active();
        for (auto &cell : dof_handler_coarse.active_cell_iterators())
          {
            if (cell->is_locally_owned())
              cell->set_active_fe_index(
                std::max((cell_other->active_fe_index() + 1u) / 2 /*bisection*/,
                         1u) -
                1);
            cell_other++;
          }
      }

      // create dof_handler
      dof_handler_fine.distribute_dofs(fe_collection);
      dof_handler_coarse.distribute_dofs(fe_collection);

      AffineConstraints<Number> constraint_coarse(
        dof_handler_coarse.locally_owned_dofs(),
        DoFTools::extract_locally_relevant_dofs(dof_handler_coarse));

      AffineConstraints<Number> constraint_fine(
        dof_handler_fine.locally_owned_dofs(),
        DoFTools::extract_locally_relevant_dofs(dof_handler_fine));

      DoFTools::make_hanging_node_constraints(dof_handler_coarse,
                                              constraint_coarse);
      constraint_coarse.close();

      DoFTools::make_hanging_node_constraints(dof_handler_fine,
                                              constraint_fine);
      constraint_fine.close();

      // setup transfer operator
      MGTwoLevelTransfer<dim, LinearAlgebra::distributed::Vector<Number>>
        transfer;
      transfer.reinit(dof_handler_fine,
                      dof_handler_coarse,
                      constraint_fine,
                      constraint_coarse);



      auto locally_relevant_dofs_coarse =
        dof_handler_coarse.locally_owned_dofs();
      locally_relevant_dofs_coarse.add_indices(
        transfer.partitioner_coarse->ghost_indices());

      DynamicSparsityPattern sparsity_restriction(dof_handler_coarse.n_dofs(),
                                                  dof_handler_fine.n_dofs(),
                                                  locally_relevant_dofs_coarse);
      transfer.make_transfer_sparsity_pattern(sparsity_restriction, false);
      sparsity_restriction
        .compress(); // does nothing here, but put here for reference as other
                     // sparsity types require this
      SparsityTools::distribute_sparsity_pattern(
        sparsity_restriction,
        transfer.partitioner_coarse->locally_owned_range(),
        dof_handler_coarse.get_communicator(),
        locally_relevant_dofs_coarse);
      TrilinosWrappers::SparseMatrix restriction_matrix;
      restriction_matrix.reinit(dof_handler_coarse.locally_owned_dofs(),
                                dof_handler_fine.locally_owned_dofs(),
                                sparsity_restriction);
      transfer.assemble_restriction_as_matrix(restriction_matrix);
      restriction_matrix.compress(VectorOperation::add);


      auto locally_relevant_dofs_fine = dof_handler_fine.locally_owned_dofs();
      locally_relevant_dofs_fine.add_indices(
        transfer.partitioner_fine->ghost_indices());

      DynamicSparsityPattern sparsity_prolongation(dof_handler_fine.n_dofs(),
                                                   dof_handler_coarse.n_dofs(),
                                                   locally_relevant_dofs_fine);
      transfer.make_transfer_sparsity_pattern(sparsity_prolongation, true);
      sparsity_prolongation
        .compress(); // does nothing here, but put here for reference as other
                     // sparsity types require this
      SparsityTools::distribute_sparsity_pattern(
        sparsity_prolongation,
        transfer.partitioner_fine->locally_owned_range(),
        dof_handler_fine.get_communicator(),
        locally_relevant_dofs_fine);
      TrilinosWrappers::SparseMatrix prolongation_matrix;
      prolongation_matrix.reinit(dof_handler_fine.locally_owned_dofs(),
                                 dof_handler_coarse.locally_owned_dofs(),
                                 sparsity_prolongation);
      transfer.assemble_prolongation_as_matrix(prolongation_matrix);
      prolongation_matrix.compress(VectorOperation::add);

      test_assembled_transfer_operator<dim, Number>(transfer,
                                                    restriction_matrix,
                                                    prolongation_matrix,
                                                    dof_handler_fine,
                                                    dof_handler_coarse);


      deallog.pop();
      // break if all cells on coarse level have active_fe_index=0
      {
        types::fe_index min = 0;

        for (auto &cell : dof_handler_coarse.active_cell_iterators())
          {
            if (cell->is_locally_owned())
              min = std::max(min, cell->active_fe_index());
          }

        min = Utilities::MPI::max(min, MPI_COMM_WORLD);

        if (min == 0)
          break;
      }
    }
}

template <int dim, typename Number>
void
test()
{
  {
    deallog.push("CG<2>");
    do_test<dim, Number>();
    deallog.pop();
  }
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    all;

  deallog.precision(8);

  test<2, double>();
}
