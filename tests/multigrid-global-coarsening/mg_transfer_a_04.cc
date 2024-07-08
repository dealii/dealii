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


/**
 * Test transfer operator for geometric coarsening without setting up MPI,
 * otherwise analogous to mg_transfer_a_01
 */

#include <deal.II/base/logstream.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/multigrid/mg_transfer_global_coarsening.h>

#include "mg_transfer_util.h"

using namespace dealii;

template <int dim, typename Number>
void
do_test(const FiniteElement<dim> &fe_fine, const FiniteElement<dim> &fe_coarse)
{
  auto create_fine_grid = [](Triangulation<dim> &tria) {
    GridGenerator::hyper_cube(tria);
    tria.refine_global();

    for (auto &cell : tria.active_cell_iterators())
      if (cell->is_active() && cell->center()[0] < 0.5)
        cell->set_refine_flag();
    tria.execute_coarsening_and_refinement();
  };

  auto execute_global_coarsening = [](Triangulation<dim> &tria) {
    for (auto &cell : tria.active_cell_iterators())
      cell->set_coarsen_flag();
    tria.execute_coarsening_and_refinement();
  };

  // create coarse grid
  Triangulation<dim> tria_coarse;
  create_fine_grid(tria_coarse);
  execute_global_coarsening(tria_coarse);

  // create fine grid
  Triangulation<dim> tria_fine;
  create_fine_grid(tria_fine);

  // setup dof-handlers
  DoFHandler<dim> dof_handler_fine(tria_fine);
  dof_handler_fine.distribute_dofs(fe_fine);

  DoFHandler<dim> dof_handler_coarse(tria_coarse);
  dof_handler_coarse.distribute_dofs(fe_coarse);

  // setup constraint matrix
  AffineConstraints<Number> constraint_coarse;
  DoFTools::make_hanging_node_constraints(dof_handler_coarse,
                                          constraint_coarse);
  constraint_coarse.close();

  AffineConstraints<Number> constraint_fine;
  DoFTools::make_hanging_node_constraints(dof_handler_fine, constraint_fine);
  constraint_fine.close();

  // setup transfer operator
  MGTwoLevelTransfer<dim, LinearAlgebra::distributed::Vector<Number>> transfer;
  transfer.reinit(dof_handler_fine,
                  dof_handler_coarse,
                  constraint_fine,
                  constraint_coarse);

  deallog << "test first time" << std::endl;
  test_transfer_operator(transfer, dof_handler_fine, dof_handler_coarse);

  deallog << "test second time" << std::endl;
  test_transfer_operator(transfer, dof_handler_fine, dof_handler_coarse);
}

template <int dim, typename Number>
void
test(int fe_degree)
{
  const auto str_fine   = std::to_string(fe_degree);
  const auto str_coarse = std::to_string(fe_degree);

  if (fe_degree > 0)
    {
      deallog.push("CG<2>(" + str_fine + ")<->CG<2>(" + str_coarse + ")");
      do_test<dim, double>(FE_Q<dim>(fe_degree), FE_Q<dim>(fe_degree));
      deallog.pop();
    }

  if (fe_degree > 0)
    {
      deallog.push("DG<2>(" + str_fine + ")<->CG<2>(" + str_coarse + ")");
      do_test<dim, double>(FE_DGQ<dim>(fe_degree), FE_Q<dim>(fe_degree));
      deallog.pop();
    }

  {
    deallog.push("DG<2>(" + str_fine + ")<->DG<2>(" + str_coarse + ")");
    do_test<dim, double>(FE_DGQ<dim>(fe_degree), FE_DGQ<dim>(fe_degree));
    deallog.pop();
  }
}

int
main()
{
  initlog();

  deallog.precision(8);

  for (unsigned int i = 0; i < 2; ++i)
    test<2, double>(i);
}
