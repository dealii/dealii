// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

// Consecutive matrix solves must replace the factorization; two-vector solves
// explicitly reuse it. Check against independently computed solutions and RHSs.
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>

#include "../tests.h"

void
set_matrix(PETScWrappers::SparseMatrix &matrix,
           const double                 a00,
           const double                 a01,
           const double                 a11)
{
  matrix.set(0, 0, a00);
  matrix.set(0, 1, a01);
  matrix.set(1, 0, a01);
  matrix.set(1, 1, a11);
  matrix.compress(VectorOperation::insert);
}

void
check(const PETScWrappers::SparseMatrix &matrix,
      const PETScWrappers::MPI::Vector  &x,
      const PETScWrappers::MPI::Vector  &b,
      const double                       x0,
      const double                       x1)
{
  AssertThrow(std::abs(x(0) - x0) < 1e-12, ExcInternalError());
  AssertThrow(std::abs(x(1) - x1) < 1e-12, ExcInternalError());
  PETScWrappers::MPI::Vector residual(x);
  AssertThrow(matrix.residual(residual, x, b) < 1e-12, ExcInternalError());
}

void
test(const bool symmetric)
{
  PETScWrappers::SparseMatrix A(2, 2, 2), B(2, 2, 2);
  set_matrix(A, 4., 1., 3.);
  set_matrix(B, 6., -1., 2.);
  IndexSet owned(2);
  owned.add_range(0, 2);
  PETScWrappers::MPI::Vector x(owned, MPI_COMM_SELF), b(owned, MPI_COMM_SELF);
  const auto                 set_rhs = [&](const double b0, const double b1) {
    b(0) = b0;
    b(1) = b1;
    b.compress(VectorOperation::insert);
  };
  SolverControl                    control;
  PETScWrappers::SparseDirectMUMPS solver(control);
  solver.set_symmetric_mode(symmetric);

  set_rhs(2., -5.);
  solver.solve(A, x, b);
  check(A, x, b, 1., -2.);

  set_rhs(8., -5.);
  solver.solve(B, x, b);
  check(B, x, b, 1., -2.);

  // Changing entries of the same matrix object must also replace its factors.
  set_matrix(B, 3., -1., 5.);
  set_rhs(5., -11.);
  solver.solve(B, x, b);
  check(B, x, b, 1., -2.);

  solver.initialize(A);
  set_rhs(2., -5.);
  solver.solve(x, b);
  check(A, x, b, 1., -2.);
  set_rhs(-10., 3.);
  solver.solve(x, b);
  check(A, x, b, -3., 2.);

  // Explicit initialization must replace the previously reusable factors too.
  solver.initialize(B);
  set_rhs(-11., 13.);
  solver.solve(x, b);
  check(B, x, b, -3., 2.);
  deallog << (symmetric ? "symmetric" : "nonsymmetric") << " OK" << std::endl;
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  initlog();
  test(false);
  test(true);
  return 0;
}
