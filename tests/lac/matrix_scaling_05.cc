// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2001 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test matrix scaling on small sparse matrix EZ

#include <deal.II/base/numbers.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/matrix_scaling.h>
#include <deal.II/lac/sparse_matrix_ez.h>
#include <deal.II/lac/vector.h>

#include "../tests.h"

using namespace dealii;

int
main(int argc, char **argv)
{
  initlog();

  unsigned int dim = 5;

  SparseMatrixEZ<double> B(dim, dim);

  for (unsigned int i = 0; i < dim; i++)
    B.set(i, i, 1.0 + i);
  ;

  MatrixScaling::AdditionalData control;
  MatrixScaling                 scaler(control);

  std::ostringstream oss1, oss2;

  deallog << "Sparse Matrix EZ: " << std::endl;

  B.print(oss1);
  deallog << oss1.str();
  scaler.find_scaling_and_scale_matrix(B);

  deallog << "Scaled Sparse Matrix EZ: " << std::endl;
  B.print(oss2);
  deallog << oss2.str();

  const Vector<double> &row_scaling    = scaler.get_row_scaling();
  const Vector<double> &column_scaling = scaler.get_column_scaling();

  deallog << "Reciprocal of scaling vectors squared " << std::endl;

  for (unsigned int i = 0; i < dim; i++)
    deallog << 1.0 / (row_scaling[i] * row_scaling[i]) << " ";
  deallog << std::endl;
  for (unsigned int i = 0; i < dim; i++)
    deallog << 1.0 / (column_scaling[i] * column_scaling[i]) << " ";
  deallog << std::endl;
}
