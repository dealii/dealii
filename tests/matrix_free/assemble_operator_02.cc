// ---------------------------------------------------------------------
//
// Copyright (C) 2014 - 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/linear_operator.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparsity_pattern.h>

#include <deal.II/matrix_free/assemble_operator.h>

#include "../tests.h"

#include "../testmatrix.h"

// Test case for MatrixFreeTools::assemble_operator.
int
main()
{
  initlog(true);

  unsigned int nx = 16;
  unsigned int ny = 16;
  unsigned int n  = nx * ny;

  FDMatrix generator(nx, ny);

  DynamicSparsityPattern dsp;
  dsp.reinit(n, n);
  generator.five_point_structure(dsp);
  dsp.compress();

  SparsityPattern pattern;
  pattern.copy_from(dsp);
  pattern.compress();

  SparseMatrix<double> M;
  M.reinit(pattern);
  generator.five_point(M);
  Assert(M.frobenius_norm() > 0, ExcInternalError());
  // Encapsulate M into a linear operator.
  // This is merely to illustrate that assemble_operator does only require
  // vmult(), as we could have used M directly.
  auto op = linear_operator(M);

  // Compute without cache.
  SparseMatrix<double> matrix;
  matrix.reinit(pattern);
  std::shared_ptr<MatrixFreeTools::GraphCache> cache;
  MatrixFreeTools::assemble_operator(matrix, op, pattern, cache);
  deallog << "n = " << n << std::endl;
  deallog << "n_vmults = " << cache->num_colors << std::endl;
  // Compute difference
  matrix.add(-1.0, M);
  deallog << std::boolalpha << (0.0 == matrix.frobenius_norm()) << std::endl;

  return 0;
}
