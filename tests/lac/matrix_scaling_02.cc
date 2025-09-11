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

// Test matrix scaling on sparse and full linear systems. Check correctness of
// solution

#include <deal.II/base/numbers.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/matrix_scaling.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/vector.h>

#include "../tests.h"

using namespace dealii;

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, numbers::invalid_unsigned_int);
  initlog();

  unsigned int           dim = 5;
  DynamicSparsityPattern dsp(dim, dim);

  for (unsigned int i = 0; i < dim; i++)
    dsp.add(i, i);

  SparsityPattern sp;
  sp.copy_from(dsp);
  SparseMatrix<double> B(sp);

  for (unsigned int i = 0; i < dim; i++)
    B.set(i, i, 1.0 + i);

  FullMatrix<double> A(dim, dim), Ainv(dim, dim);
  A.copy_from(B);
  Ainv.invert(A);

  Vector<double> x(dim), rhs(dim);

  for (unsigned int i = 0; i < dim; i++)
    rhs[i] = 1.0;

  MatrixScaling::AdditionalData control;
  MatrixScaling                 scaler(control);

  deallog << "Sparse Matrix: " << std::endl;

  SparseDirectUMFPACK Binv;
  Binv.initialize(B);
  Binv.vmult(x, rhs);
  deallog << "Solution of Ax=b: " << std::endl;
  for (unsigned int i = 0; i < dim; i++)
    deallog << x[i] << " ";
  deallog << std::endl;

  scaler.find_scaling_and_scale_linear_system(B, rhs);
  scaler.scale_system_solution(rhs);

  deallog << "Solution of Ax=b scaled: " << std::endl;
  for (unsigned int i = 0; i < dim; i++)
    deallog << rhs[i] << " ";
  deallog << std::endl;

  deallog << "Full Matrix: " << std::endl;

  for (unsigned int i = 0; i < dim; i++)
    rhs[i] = 1.0;

  Ainv.vmult(x, rhs);
  deallog << "Solution of Ax=b: " << std::endl;
  for (unsigned int i = 0; i < dim; i++)
    deallog << x[i] << " ";
  deallog << std::endl;

  scaler.find_scaling_and_scale_linear_system(A, rhs);
  scaler.scale_system_solution(rhs);

  deallog << "Solution of Ax=b scaled: " << std::endl;
  for (unsigned int i = 0; i < dim; i++)
    deallog << rhs[i] << " ";
  deallog << std::endl;
}
