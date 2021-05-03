// ---------------------------------------------------------------------
//
// Copyright (C) 2020 by the deal.II authors
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


/**
 * Like mg_transfer_a_01.cc but working with FECollections with more than one
 * FE.
 */

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>

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
  const auto create_fine_grid = [](Triangulation<dim> &tria) {
    GridGenerator::hyper_cube(tria);
    tria.refine_global();

    for (auto &cell : tria.active_cell_iterators())
      if (cell->active() && cell->center()[0] < 0.5)
        cell->set_refine_flag();
    tria.execute_coarsening_and_refinement();
  };

  const auto execute_global_coarsening = [](Triangulation<dim> &tria) {
    for (auto &cell : tria.active_cell_iterators())
      cell->set_coarsen_flag();
    tria.execute_coarsening_and_refinement();
  };

  const hp::FECollection<dim> fe(fe_fine, fe_coarse);

  // create coarse grid
  parallel::distributed::Triangulation<dim> tria_coarse(MPI_COMM_WORLD);
  create_fine_grid(tria_coarse);
  execute_global_coarsening(tria_coarse);

  // create fine grid
  parallel::distributed::Triangulation<dim> tria_fine(MPI_COMM_WORLD);
  create_fine_grid(tria_fine);

  // setup dof-handlers
  DoFHandler<dim> dof_handler_fine(tria_fine);
  for (const auto &cell : dof_handler_fine.active_cell_iterators())
    if (cell->is_locally_owned())
      cell->set_active_fe_index(0);
  dof_handler_fine.distribute_dofs(fe);

  DoFHandler<dim> dof_handler_coarse(tria_coarse);
  for (const auto &cell : dof_handler_coarse.active_cell_iterators())
    if (cell->is_locally_owned())
      cell->set_active_fe_index(1);
  dof_handler_coarse.distribute_dofs(fe);

  // setup constraint matrix
  AffineConstraints<Number> constraint_coarse;
  DoFTools::make_hanging_node_constraints(dof_handler_coarse,
                                          constraint_coarse);
  constraint_coarse.close();

  AffineConstraints<Number> constraint_fine;
  DoFTools::make_hanging_node_constraints(dof_handler_fine, constraint_fine);
  constraint_coarse.close();

  // setup transfer operator
  MGTwoLevelTransfer<dim, LinearAlgebra::distributed::Vector<Number>> transfer;
  transfer.reinit_geometric_transfer(dof_handler_fine,
                                     dof_handler_coarse,
                                     constraint_fine,
                                     constraint_coarse);

  test_transfer_operator(transfer, dof_handler_fine, dof_handler_coarse);
}

template <int dim, typename Number>
void
test(int fe_degree)
{
  const auto str_fine   = std::to_string(fe_degree);
  const auto str_coarse = std::to_string(fe_degree);

  {
    deallog.push("CG<2>(" + str_fine + ")<->CG<2>(" + str_coarse + ")");
    do_test<dim, double>(FE_Q<dim>(fe_degree), FE_Q<dim>(fe_degree));
    deallog.pop();
  }

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
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    all;

  deallog.precision(8);

  for (unsigned int i = 1; i < 5; i++)
    test<2, double>(i);
}
