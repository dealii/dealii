// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// tests that GMRES builds an orthonormal basis properly for a larger matrix
// and basis size than gmres_orthogonalize_01

#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include "../tests.h"



template <typename number>
void
test()
{
  const unsigned int n = 200;
  Vector<number>     rhs(n), sol(n);
  rhs = 1.;

  // only add diagonal entries
  SparsityPattern sp(n, n);
  sp.compress();
  SparseMatrix<number> matrix(sp);

  for (unsigned int i = 0; i < n; ++i)
    matrix.diag_element(i) = (i + 1);

  SolverControl                                        control(1000,
                        1e3 * std::numeric_limits<number>::epsilon(),
                        false,
                        true);
  typename SolverGMRES<Vector<number>>::AdditionalData data;
  data.max_basis_size = 200;
  data.orthogonalization_strategy =
    LinearAlgebra::OrthogonalizationStrategy::modified_gram_schmidt;

  SolverGMRES<Vector<number>> solver(control, data);
  auto print_re_orthogonalization = [](int accumulated_iterations) {
    deallog.get_file_stream() << "Re-orthogonalization enabled at step "
                              << accumulated_iterations << std::endl;
  };
  solver.connect_re_orthogonalization_slot(print_re_orthogonalization);
  solver.solve(matrix, sol, rhs, PreconditionIdentity());
}

int
main()
{
  initlog();
  deallog << std::setprecision(4);

  deallog.push("double");
  test<double>();
  deallog.pop();
  deallog.push("float");
  test<float>();
  deallog.pop();
}
