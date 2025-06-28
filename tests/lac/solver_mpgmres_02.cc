// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Check that we can call FGMRES::solve() and MPGMRES::solve() without a
// preconditioner

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/sparse_matrix.h>

#include "../tests.h"

int
main()
{
  initlog();

  const unsigned int N = 500;
  const unsigned int M = 2 * N;

  SparsityPattern sparsity_pattern(M, M, 2);
  for (unsigned int i = 0; i < M - 1; ++i)
    sparsity_pattern.add(i, i + 1);
  sparsity_pattern.compress();
  SparseMatrix<double> matrix(sparsity_pattern);

  for (unsigned int i = 0; i < N; ++i)
    matrix.diag_element(i) = 1.0;
  for (unsigned int i = N; i < M; ++i)
    matrix.diag_element(i) = 1.;

  Vector<double> rhs(M);
  Vector<double> sol(M);
  rhs = 1.;

  SolverControl control(M, 1.e-2);

  deallog.get_file_stream() << "FGMRES:\n";
  SolverFGMRES<Vector<double>> solver_fgmres(control);
  check_solver_within_range(solver_fgmres.solve(matrix, sol, rhs),
                            control.last_step(),
                            1,
                            10);

  deallog.get_file_stream() << "\nMPGMRES:\n";
  sol = 0.;

  SolverMPGMRES<Vector<double>> solver_mpgmres(control);
  check_solver_within_range(solver_mpgmres.solve(matrix, sol, rhs),
                            control.last_step(),
                            1,
                            10);
}
