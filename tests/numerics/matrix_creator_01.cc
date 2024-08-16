// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test MatrixCreator::create_mass_matrix() and
// MatrixCreator::create_laplace_matrix() for Trilinos and PETSc matrices in
// parallel.

#include <deal.II/base/mpi.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/generic_linear_algebra.h>
#include <deal.II/lac/sparsity_tools.h>

#include <deal.II/numerics/matrix_tools.h>

#include "../tests.h"


template <int dim, typename MatrixType>
void
test()
{
  parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(tria);
  tria.refine_global(2);

  FE_Q<dim>   fe(2);
  QGauss<dim> quad(3);

  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  const auto      mpi_comm   = dof_handler.get_mpi_communicator();
  const IndexSet &owned_dofs = dof_handler.locally_owned_dofs();

  const IndexSet relevant_dofs =
    DoFTools::extract_locally_relevant_dofs(dof_handler);

  AffineConstraints<typename MatrixType::value_type> constraints;

  DynamicSparsityPattern sparsity_pattern(relevant_dofs);
  DoFTools::make_sparsity_pattern(dof_handler, sparsity_pattern, constraints);
  SparsityTools::distribute_sparsity_pattern(sparsity_pattern,
                                             owned_dofs,
                                             mpi_comm,
                                             relevant_dofs);

  MatrixType matrix;
  matrix.reinit(owned_dofs, owned_dofs, sparsity_pattern, mpi_comm);

  MatrixCreator::create_laplace_matrix(dof_handler, quad, matrix);
  deallog << "Laplace matrix: " << matrix.frobenius_norm() << std::endl;

  matrix = 0.0;
  MatrixCreator::create_mass_matrix(dof_handler, quad, matrix);
  deallog << "Mass matrix:    " << matrix.frobenius_norm() << std::endl;
}

int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  mpi_initlog();

  {
    deallog.push("Trilinos");
    test<2, TrilinosWrappers::SparseMatrix>();
    deallog.pop();
  }
  {
    deallog.push("PETSc");
    test<2, PETScWrappers::MPI::SparseMatrix>();
    deallog.pop();
  }
}
