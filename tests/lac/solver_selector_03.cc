// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test the SolverSelector class using the constructor that initializes the
// object already.

#include <deal.II/lac/solver_selector.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/vector_memory.h>

#include "../tests.h"

#include "../testmatrix.h"


template <typename MatrixType, typename VectorType>
void
check(const MatrixType &A, const VectorType &f)
{
  std::vector<std::string> names;
  names.push_back("cg");
  names.push_back("gmres");

  SolverControl                          cont(100, 1.e-7, false, true);
  SolverSelector<VectorType>             solver;
  PreconditionSSOR<SparseMatrix<double>> pre;
  pre.initialize(A);

  VectorType u;
  u.reinit(f);

  for (const auto &name : names)
    {
      SolverSelector<VectorType> solver(name, cont);
      u = 0.;
      solver.solve(A, u, f, pre);
    }
}


int
main()
{
  initlog();
  deallog << std::setprecision(4);

  unsigned int size = 37;
  unsigned int dim  = (size - 1) * (size - 1);

  deallog << "Size " << size << " Unknowns " << dim << std::endl;

  // Make matrix
  FDMatrix        testproblem(size, size);
  SparsityPattern structure(dim, dim, 5);
  testproblem.five_point_structure(structure);
  structure.compress();
  SparseMatrix<double> A(structure);
  testproblem.five_point(A);
  Vector<double> f(dim);
  f = 1.;

  check(A, f);
}
