// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2012 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// test lucky breakdown in GMRES (and others)

#include <deal.II/base/point.h>

#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/solver_idr.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/vector_memory.h>

#include "../tests.h"

#include "../testmatrix.h"

template <typename SolverType>
void
test()
{
  const unsigned int size = 3;
  SparsityPattern    sparsity(size, size, 1);
  sparsity.compress();
  SparseMatrix<double> mat;
  mat.reinit(sparsity);
  mat = IdentityMatrix(size);

  Vector<double> rhs;
  Vector<double> solvec;
  solvec.reinit(size);

  rhs.reinit(size);
  rhs(size - 1) = 1.0;

  SolverControl solvctrl(1000, 1e-12, true);
  SolverType    solver(solvctrl);

  PreconditionIdentity precond;
  solver.solve(mat, solvec, rhs, precond);
  solvec.print(deallog.get_file_stream());
}

int
main()
{
  initlog();
  deallog << std::setprecision(4);

  test<SolverGMRES<Vector<double>>>();
  test<SolverCG<Vector<double>>>();
  test<SolverFGMRES<Vector<double>>>();
  test<SolverIDR<Vector<double>>>();
}
