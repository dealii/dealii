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
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/linear_operator.h>

#include <deal.II/matrix_free/assemble_operator.h>

#include "../tests.h"

// Small test case for MatrixFreeTools::assemble_operator.
// Uses a FullMatrix as operator as well as for the result.
// The pattern of the result is not important,
// as long as it is a superset of the provided sparsity pattern.
int
main()
{
  initlog(true);

  DynamicSparsityPattern pattern;
  pattern.reinit(8, 6);
  pattern.add(0, 0);
  pattern.add(0, 1);
  pattern.add(0, 5);
  pattern.add(1, 0);
  pattern.add(1, 1);
  pattern.add(1, 2);
  pattern.add(2, 1);
  pattern.add(2, 2);
  pattern.add(2, 3);
  pattern.add(2, 5);
  pattern.add(3, 2);
  pattern.add(3, 3);
  pattern.add(4, 4);
  pattern.add(5, 0);
  pattern.add(5, 2);
  pattern.add(5, 5);
  pattern.add(6, 4);
  pattern.add(6, 5);
  pattern.add(7, 5);
  pattern.compress();

  FullMatrix<double> M;
  M.reinit(8, 6);
  for (auto p : pattern)
    M(p.row(), p.column()) = std::sin(123.0 + p.row() * p.column());
  //   M.print_formatted(std::cout, 2, false, 4);

  // Encapsulate M into a linear operator.
  // This is merely to illustrate that assemble_operator does only require
  // vmult(), as we could have used M directly.
  auto op = linear_operator(M);

  // Compute without cache.
  FullMatrix<double> matrix;
  matrix.reinit(8, 6);
  std::shared_ptr<MatrixFreeTools::GraphCache> cache;
  MatrixFreeTools::assemble_operator(matrix, op, pattern, cache);
  deallog << "n = " << matrix.n_rows() << std::endl;
  deallog << "n_vmults = " << cache->num_colors << std::endl;
  //   matrix.print_formatted(deallog, 2, false, 4);
  deallog << std::boolalpha << (M == matrix) << std::endl;

  // Reuse cache.
  FullMatrix<double> matrix2;
  matrix2.reinit(8, 6);
  MatrixFreeTools::assemble_operator(matrix2, op, pattern, cache);
  //   matrix2.print_formatted(deallog, 2, false, 4);
  deallog << std::boolalpha << (M == matrix2) << std::endl;

  return 0;
}
