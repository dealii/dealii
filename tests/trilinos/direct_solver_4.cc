// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// tests Trilinos direct solvers for multi vectors

#include <deal.II/base/function.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_solver.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"

int
main(int argc, char **argv)
{
  initlog();

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv);

  const unsigned int n = 5;
  const unsigned int m = 3;

  // initialize the matrix
  FullMatrix<double> A(n, n);

  for (unsigned int i = 0; i < n; ++i)
    A[i][i] = 2.0;

  for (unsigned int i = 0; i < n - 1; ++i)
    {
      A[i][i + 1] = -1.0;
      A[i + 1][i] = -1.0;
    }

  // initialize the right-hand-side vector
  FullMatrix<double> B(n, m);

  for (unsigned int i = 0; i < m; ++i)
    B[i][i] = 1;

  // solve with direct solver
  {
    auto A_copy = A;
    A_copy.gauss_jordan();

    FullMatrix<double> X(n, m);
    A_copy.mmult(X, B);

    X.print(deallog.get_file_stream());
  }

  deallog << std::endl;

  // solve with sparse direct solver
  {
    std::vector<std::pair<types::global_dof_index, types::global_dof_index>>
      entries;

    for (unsigned int i = 0; i < A.m(); ++i)
      for (unsigned int j = 0; j < A.n(); ++j)
        if (A[i][j] != 0)
          entries.emplace_back(i, j);

    SparsityPattern sp(n, n);
    sp.add_entries(entries);
    sp.compress();

    TrilinosWrappers::SparseMatrix A_sparse;
    A_sparse.reinit(sp);

    for (unsigned int i = 0; i < A.m(); ++i)
      for (unsigned int j = 0; j < A.n(); ++j)
        if (A[i][j] != 0)
          A_sparse.set(i, j, A[i][j]);

    TrilinosWrappers::SolverDirect solver;

    FullMatrix<double> X(n, m);
    solver.solve(A_sparse, X, B);
    X.print(deallog.get_file_stream());
  }
}
