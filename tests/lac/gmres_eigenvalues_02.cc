// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2018 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


// test that GMRES eigenvalue approximation doesn't crash
// when no iterations are performed.

#include <deal.II/lac/lapack_full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/vector.h>

#include "../tests.h"

template <typename number>
void
test()
{
  const unsigned int n = 2;
  Vector<number>     rhs(n), sol(n);
  rhs = 0.;

  LAPACKFullMatrix<number> matrix(n, n);

  for (unsigned int i = 0; i < n; ++i)
    matrix(i, i) = std::sqrt(i + 1);

  SolverControl control;
  {
    SolverGMRES<Vector<number>> solver(control);
    solver.connect_eigenvalues_slot(
      [](const std::vector<std::complex<double>> &eigenvalues) {
        deallog << "n_eigenvalues: " << eigenvalues.size();
      });
    solver.solve(matrix, sol, rhs, PreconditionIdentity());
  }

  {
    SolverCG<Vector<number>> solver(control);
    solver.connect_eigenvalues_slot([](const std::vector<double> &eigenvalues) {
      deallog << "n_eigenvalues: " << eigenvalues.size();
    });
    solver.solve(matrix, sol, rhs, PreconditionIdentity());
  }

  deallog << "OK" << std::endl;
}

int
main()
{
  initlog();
  test<double>();
}
