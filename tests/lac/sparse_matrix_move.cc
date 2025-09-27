// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



#include <deal.II/lac/sparse_matrix.h>

#include "../tests.h"

#include "../testmatrix.h"


void
graph_laplacian(const SparsityPattern &sparsity, SparseMatrix<double> &matrix)
{
  matrix = 0.0;

  for (SparsityPattern::const_iterator it = sparsity.begin();
       it != sparsity.end();
       ++it)
    {
      const auto i = (*it).row();
      const auto j = (*it).column();
      matrix.add(i, j, -1);
      matrix.add(i, i, 1);
    }
}


SparseMatrix<double>
graph_laplacian(const SparsityPattern &sparsity)
{
  SparseMatrix<double> A(sparsity);
  graph_laplacian(sparsity, A);

  return A;
}



int
main()
{
  initlog();
  deallog << std::setprecision(3);

  const unsigned int size = 5;

  FDMatrix     testproblem(size, size);
  unsigned int dim = (size - 1) * (size - 1);

  SparsityPattern sparsity(dim, dim, size);
  testproblem.five_point_structure(sparsity);
  sparsity.compress();

  Vector<double> x(dim), y(dim);

  {
    // Return a sparse matrix, possibly using RVO or move constructor
    SparseMatrix<double> A = graph_laplacian(sparsity);
    deallog << A.n_nonzero_elements() << std::endl;

    x = 1.0;
    y = 1.0;
    A.vmult(y, x);

    deallog << y.l2_norm() << std::endl;
  }

  {
    // Return a sparse matrix using the move constructoir
    SparseMatrix<double> B = graph_laplacian(sparsity);
    SparseMatrix<double> A = std::move(B);
    deallog << A.n_nonzero_elements() << std::endl;
    deallog << B.empty() << std::endl;
    y = 1.0;
    A.vmult(y, x);

    deallog << y.l2_norm() << std::endl;

    // Explicitly move a sparse matrix
    B = std::move(A);
    deallog << B.m() << std::endl;
    deallog << A.empty() << std::endl;
  }

  return 0;
}
