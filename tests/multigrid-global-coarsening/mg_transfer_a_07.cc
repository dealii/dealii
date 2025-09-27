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
 * Test MGTwoLevelTransfer for FE_DGP.
 */

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgp.h>
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

template <int dim, typename Number, typename MeshType>
void
test_interpolation(
  const MGTwoLevelTransferBase<LinearAlgebra::distributed::Vector<Number>>
                              &transfer,
  const MeshType              &dof_handler_fine,
  const MeshType              &dof_handler_coarse,
  const Function<dim, Number> &function,
  const unsigned int           mg_level_fine   = numbers::invalid_unsigned_int,
  const unsigned int           mg_level_coarse = numbers::invalid_unsigned_int)
{
  AffineConstraints<Number> constraint_fine(
    dof_handler_fine.locally_owned_dofs(),
    DoFTools::extract_locally_relevant_dofs(dof_handler_fine));
  DoFTools::make_hanging_node_constraints(dof_handler_fine, constraint_fine);
  constraint_fine.close();

  // perform interpolation
  LinearAlgebra::distributed::Vector<Number> src, dst;

  initialize_dof_vector(dst, dof_handler_fine, mg_level_fine);
  initialize_dof_vector(src, dof_handler_coarse, mg_level_coarse);

  // test interpolate
  AffineConstraints<Number> dummy;
  dummy.close();
  VectorTools::project(dof_handler_fine, dummy, QGauss<dim>(4), function, dst);

  transfer.interpolate(src, dst);
  print_if_non_zero(dst, 1e-7);
  print_if_non_zero(src, 1e-7);

  // prolongate back
  dst = 0.0;
  transfer.prolongate_and_add(dst, src);
  print_if_non_zero(dst, 1e-7);
  deallog << std::endl;
}

template <int dim, typename Number>
void
do_test(const FiniteElement<dim>    &fe_fine,
        const FiniteElement<dim>    &fe_coarse,
        const Function<dim, Number> &function,
        const bool                   refine_fine_mesh)
{
  // create coarse grid
  Triangulation<dim> tria_coarse;
  GridGenerator::hyper_cube(tria_coarse);

  // create fine grid
  Triangulation<dim> tria_fine;
  GridGenerator::hyper_cube(tria_fine);

  if (refine_fine_mesh)
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

  test_interpolation(transfer, dof_handler_fine, dof_handler_coarse, function);
}

template <int dim, typename Number>
void
test(int fe_degree, const Function<dim, Number> &function)
{
  const auto str_fine   = std::to_string(fe_degree);
  const auto str_coarse = std::to_string(fe_degree);
  const auto str_dim    = std::to_string(dim);

  deallog.push("FE_DGP<" + str_dim + ">(" + str_fine + ")");
  do_test<dim, double>(FE_DGP<dim>(fe_degree),
                       FE_DGP<dim>(fe_degree),
                       function,
                       true);
  do_test<dim, double>(FE_DGP<dim>(fe_degree + 1),
                       FE_DGP<dim>(fe_degree),
                       function,
                       false);
  deallog.pop();
}

template <int dim>
void
test()
{
  for (unsigned int i = 0; i <= 3; ++i)
    test<dim, double>(i, Functions::ConstantFunction<dim, double>(1.));

  deallog << std::endl << std::endl << std::endl;

  for (unsigned int i = 0; i <= 3; ++i)
    test<dim, double>(i, RightHandSideFunction<dim>());
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    all;

  deallog.precision(8);

  test<1>();
  test<2>();
}
