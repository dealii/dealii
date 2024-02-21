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
 * Test transfer operator for polynomial coarsening.
 *
 * Example:
 *
 * +--+--+------+      +--+--+------+
 * |kf|kf|      |      |kc|kc|      |
 * |--+--|  kf  |      |--+--|  kc  |
 * |kf|kf|      |  0   |kc|kc|      |
 * +--+--+------+  ->  +--+--+------+
 * |kf|kf|      |      |kc|kc|      |
 * |--+--|  kf  |      |--+--|  kc  |
 * |kf|kf|      |      |kc|kc|      |
 * +--+--+------+      +--+--+------+
 *
 *                 ... with fe_degree in the cells
 */

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/mpi.h>
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

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <deal.II/multigrid/mg_constrained_dofs.h>
#include <deal.II/multigrid/mg_transfer_global_coarsening.h>

#include "mg_transfer_util.h"

using namespace dealii;

template <int dim, typename Number>
void
do_test(const FiniteElement<dim> &fe_fine, const FiniteElement<dim> &fe_coarse)
{
  parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);

  // create grid
  GridGenerator::hyper_cube(tria);
  tria.refine_global();

  for (auto &cell : tria.active_cell_iterators())
    if (cell->is_active() && cell->center()[0] < 0.5)
      cell->set_refine_flag();
  tria.execute_coarsening_and_refinement();

  // setup dof-handlers
  DoFHandler<dim> dof_handler_fine(tria);
  dof_handler_fine.distribute_dofs(fe_fine);

  DoFHandler<dim> dof_handler_coarse(tria);
  dof_handler_coarse.distribute_dofs(fe_coarse);

  AffineConstraints<Number> constraint_coarse(
    dof_handler_coarse.locally_owned_dofs(),
    DoFTools::extract_locally_relevant_dofs(dof_handler_coarse));

  AffineConstraints<Number> constraint_fine(
    dof_handler_fine.locally_owned_dofs(),
    DoFTools::extract_locally_relevant_dofs(dof_handler_fine));

  DoFTools::make_hanging_node_constraints(dof_handler_coarse,
                                          constraint_coarse);
  constraint_coarse.close();

  DoFTools::make_hanging_node_constraints(dof_handler_fine, constraint_fine);
  constraint_fine.close();

  // setup transfer operator
  MGTwoLevelTransfer<dim, LinearAlgebra::distributed::Vector<Number>> transfer;
  transfer.reinit(dof_handler_fine,
                  dof_handler_coarse,
                  constraint_fine,
                  constraint_coarse);

  test_transfer_operator(transfer, dof_handler_fine, dof_handler_coarse);
}

template <int dim, typename Number>
void
test(int fe_degree_fine, int fe_degree_coarse)
{
  const auto str_fine   = std::to_string(fe_degree_fine);
  const auto str_coarse = std::to_string(fe_degree_coarse);

  {
    deallog.push("CG<2>(" + str_fine + ")<->CG<2>(" + str_coarse + ")");
    do_test<dim, Number>(FE_Q<dim>(fe_degree_fine),
                         FE_Q<dim>(fe_degree_coarse));
    deallog.pop();
  }

  {
    deallog.push("DG<2>(" + str_fine + ")<->CG<2>(" + str_coarse + ")");
    do_test<dim, Number>(FE_DGQ<dim>(fe_degree_fine),
                         FE_Q<dim>(fe_degree_coarse));
    deallog.pop();
  }

  {
    deallog.push("DG<2>(" + str_fine + ")<->DG<2>(" + str_coarse + ")");
    do_test<dim, Number>(FE_DGQ<dim>(fe_degree_fine),
                         FE_DGQ<dim>(fe_degree_coarse));
    deallog.pop();
  }
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    all;

  deallog.precision(8);

  for (unsigned int fe_degree_fine = 1; fe_degree_fine <= 5; ++fe_degree_fine)
    for (unsigned int fe_degree_coarse = 1; fe_degree_coarse <= fe_degree_fine;
         fe_degree_coarse++)
      test<2, double>(fe_degree_fine, fe_degree_coarse);
}
