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

// Test matrix scaling with SK algorithm on small sparse and full matrices.
// Convergence is checked and scaling vectors are printed.

#include <deal.II/base/numbers.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/matrix_scaling.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/vector.h>

#include "../tests.h"

using namespace dealii;

int
main(int argc, char **argv)
{
  initlog();

  bool converged;

  unsigned int           dim = 5;
  DynamicSparsityPattern dsp(dim, dim);

  for (unsigned int i = 0; i < dim; i++)
    dsp.add(i, i);

  SparsityPattern sp;
  sp.copy_from(dsp);
  SparseMatrix<double> B(sp);

  for (unsigned int i = 0; i < dim; i++)
    B.set(i, i, 1.0 + i);

  FullMatrix<double> A(dim, dim);
  A.copy_from(B);

  MatrixScaling::AdditionalData control;

  control.algorithm =
    MatrixScaling::AdditionalData::ScalingAlgorithm::sinkhorn_knopp;
  MatrixScaling scaler(control);

  deallog << "Sparse Matrix: " << std::endl;

  B.print(deallog);
  converged = scaler.find_scaling_and_scale_matrix(B);

  deallog << "Converged? " << (converged ? "True" : "False") << std::endl;


  deallog << "Scaled Sparse Matrix: " << std::endl;
  B.print(deallog);

  const Vector<double> &row_scaling    = scaler.get_row_scaling();
  const Vector<double> &column_scaling = scaler.get_column_scaling();

  deallog << "Reciprocal of scaling vectors " << std::endl;

  for (unsigned int i = 0; i < dim; i++)
    deallog << 1.0 / (row_scaling[i]) << " ";
  deallog << std::endl;
  for (unsigned int i = 0; i < dim; i++)
    deallog << 1.0 / (column_scaling[i]) << " ";
  deallog << std::endl;

  deallog << "Full Matrix: " << std::endl;

  A.print(deallog);
  converged = scaler.find_scaling_and_scale_matrix(A);

  deallog << "Converged? " << (converged ? "True" : "False") << std::endl;

  deallog << "Scaled Full Matrix: " << std::endl;
  A.print(deallog);

  deallog << "Reciprocal of scaling vectors " << std::endl;

  for (unsigned int i = 0; i < dim; i++)
    deallog << 1.0 / (row_scaling[i]) << " ";
  deallog << std::endl;
  for (unsigned int i = 0; i < dim; i++)
    deallog << 1.0 / (column_scaling[i]) << " ";
  deallog << std::endl;
}
