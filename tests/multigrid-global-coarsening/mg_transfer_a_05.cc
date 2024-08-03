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
 * Test prolongation with uniformly refined meshes.
 *
 * To compare with non_nested_transfer_04.
 */

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_simplex_p.h>
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

template <int dim>
class RightHandSideFunction : public Function<dim>
{
public:
  RightHandSideFunction()
    : Function<dim>(1)
  {}

  virtual double
  value(const Point<dim> &p, const unsigned int component = 0) const
  {
    (void)component;
    return p[0];
  }
};

template <int dim, typename Number>
void
do_test(const FiniteElement<dim>    &fe_fine,
        const FiniteElement<dim>    &fe_coarse,
        const Function<dim, Number> &function)
{
  // create coarse grid
  Triangulation<dim> tria_coarse; // TODO
  if (fe_coarse.reference_cell().is_simplex())
    GridGenerator::subdivided_hyper_cube_with_simplices(tria_coarse, 1);
  else
    GridGenerator::hyper_cube(tria_coarse);
  tria_coarse.refine_global();

  // create fine grid
  Triangulation<dim> tria_fine; // TODO

  if (fe_fine.reference_cell().is_simplex())
    GridGenerator::subdivided_hyper_cube_with_simplices(tria_fine, 1);
  else
    GridGenerator::hyper_cube(tria_fine);
  tria_fine.refine_global();
  tria_fine.refine_global();

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
  MGTwoLevelTransfer<dim, LinearAlgebra::distributed::Vector<Number>> transfer;

  transfer.reinit(dof_handler_fine,
                  dof_handler_coarse,
                  constraint_fine,
                  constraint_coarse);

  test_non_nested_transfer(transfer,
                           dof_handler_fine,
                           dof_handler_coarse,
                           function);
}

template <int dim, typename Number>
void
test(int                          fe_degree,
     const Function<dim, Number> &function,
     const bool                   do_simplex_mesh)
{
  const auto str_fine   = std::to_string(fe_degree);
  const auto str_coarse = std::to_string(fe_degree);
  const auto str_dim    = std::to_string(dim);

  if ((fe_degree > 0) && (do_simplex_mesh == false))
    {
      deallog.push("HCG<" + str_dim + ">(" + str_fine + ")<->HCG<" + str_dim +
                   ">(" + str_coarse + ")");
      do_test<dim, double>(FE_Q<dim>(fe_degree),
                           FE_Q<dim>(fe_degree),
                           function);
      deallog.pop();
    }

  if ((fe_degree > 0) && do_simplex_mesh)
    {
      deallog.push("SCG<" + str_dim + ">(" + str_fine + ")<->SCG<" + str_dim +
                   ">(" + str_coarse + ")");
      do_test<dim, double>(FE_SimplexP<dim>(fe_degree),
                           FE_SimplexP<dim>(fe_degree),
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

  // Functions::ConstantFunction<2, double> fu(1.);
  RightHandSideFunction<2> fu;
  RightHandSideFunction<3> fu3;

  for (unsigned int i = 0; i <= 4; ++i)
    test<2, double>(i, fu, false);

  for (unsigned int i = 0; i <= 3; ++i)
    test<2, double>(i, fu, true);

  for (unsigned int i = 1; i <= 3; ++i)
    test<3, double>(i, fu3, true);
}
