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
 * Test prolongation and restriction with "non-nested" triangulations.
 *
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
#include <deal.II/grid/grid_tools.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <deal.II/multigrid/mg_constrained_dofs.h>
#include <deal.II/multigrid/mg_transfer_global_coarsening.h>

#include "mg_transfer_util.h"

using namespace dealii;

template <int dim, typename Number>
void
do_test(const FiniteElement<dim>    &fe_fine,
        const FiniteElement<dim>    &fe_coarse,
        const Function<dim, Number> &function)
{
  // create coarse grid
  parallel::distributed::Triangulation<dim> tria_coarse(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(tria_coarse, -1. + .2, 1. + .2);
  const double angle = .5 * numbers::PI_4;
  GridTools::rotate(angle, tria_coarse);
  tria_coarse.refine_global(2);

  // create fine grid
  parallel::distributed::Triangulation<dim> tria_fine(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(tria_fine, -2., +2.);
  tria_fine.refine_global(3);

  // setup dof-handlers
  DoFHandler<dim> dof_handler_fine(tria_fine);
  dof_handler_fine.distribute_dofs(fe_fine);

  DoFHandler<dim> dof_handler_coarse(tria_coarse);
  dof_handler_coarse.distribute_dofs(fe_coarse);

  // setup constraint matrix
  AffineConstraints<Number> constraint_coarse;
  constraint_coarse.close();

  AffineConstraints<Number> constraint_fine;
  constraint_fine.close();

  // setup transfer operator
  typename MGTwoLevelTransferNonNested<
    dim,
    LinearAlgebra::distributed::Vector<Number>>::AdditionalData data;
  data.enforce_all_points_found = false;
  MGTwoLevelTransferNonNested<dim, LinearAlgebra::distributed::Vector<Number>>
                 transfer(data);
  MappingQ1<dim> mapping_fine, mapping_coarse;
  transfer.reinit(dof_handler_fine,
                  dof_handler_coarse,
                  mapping_fine,
                  mapping_coarse,
                  constraint_fine,
                  constraint_coarse);

  test_non_nested_transfer(transfer,
                           dof_handler_fine,
                           dof_handler_coarse,
                           function);
}

template <int dim, typename Number>
void
test(int fe_degree, const Function<dim, Number> &function)
{
  const auto str_fine   = std::to_string(fe_degree);
  const auto str_coarse = std::to_string(fe_degree);

  if (fe_degree > 0)
    {
      deallog.push("CG<2>(" + str_fine + ")<->CG<2>(" + str_coarse + ")");
      do_test<dim, double>(FE_Q<dim>(fe_degree),
                           FE_Q<dim>(fe_degree),
                           function);
      deallog.pop();
    }

  if (fe_degree > 0)
    {
      deallog.push("DG<2>(" + str_fine + ")<->CG<2>(" + str_coarse + ")");
      do_test<dim, double>(FE_DGQ<dim>(fe_degree),
                           FE_Q<dim>(fe_degree),
                           function);
      deallog.pop();
    }

  {
    deallog.push("DG<2>(" + str_fine + ")<->DG<2>(" + str_coarse + ")");
    do_test<dim, double>(FE_DGQ<dim>(fe_degree),
                         FE_DGQ<dim>(fe_degree),
                         function);
    deallog.pop();
  }
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    all;

  deallog.precision(8);
  // Functions::Monomial<2> linear_function(Tensor<1, 2>({1, 0}));   // f(x,y)=
  // x
  Functions::Monomial<2> bilinear_function(Tensor<1, 2>({1, 1})); // f(x,y)= xy

  for (unsigned int i = 0; i < 3; ++i)
    test<2, double>(i, bilinear_function);
}
