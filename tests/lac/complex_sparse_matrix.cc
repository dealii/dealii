// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test SparseMatrix<complex<number>> copy_from, add from SparseMatrix<number>,
// and vmult, vmult_add, ... for mixed real and complex vectors

#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/vector.h>

#include <complex>

#include "../tests.h"


template <typename number>
void
check_operations()
{
  using namespace std::complex_literals;

  SparsityPattern pattern(1, 1, 1);
  pattern.add(0, 0);
  pattern.compress();
  SparseMatrix<number>               real_matrix(pattern);
  SparseMatrix<std::complex<number>> cmpl_matrix(pattern);

  real_matrix.set(0, 0, 1.0);
  deallog << "real_matrix.set(0, 0, 1.0) "
          << (real_matrix(0, 0) == 1.0 ? "True" : "False") << std::endl;


  cmpl_matrix.copy_from(real_matrix);
  Assert((cmpl_matrix(0, 0).real() == real_matrix(0, 0)) &&
           (cmpl_matrix(0, 0).imag() == 0.0),
         ExcInternalError());
  deallog << "cmpl_matrix.copy_from(real_matrix) True" << std::endl;


  cmpl_matrix.add(0, 0, 1i);
  Assert((cmpl_matrix(0, 0).real() == 1.0) && (cmpl_matrix(0, 0).imag() == 1.0),
         ExcInternalError());
  deallog << "cmpl_matrix.add(0, 0, 1i) True" << std::endl;

  Vector<number>               real_vec({2.});
  Vector<std::complex<number>> cmpl_vec({0.});

  cmpl_matrix.vmult(cmpl_vec, real_vec);
  Assert((cmpl_vec(0).real() == 2.0) && (cmpl_vec(0).imag() == 2.0),
         ExcInternalError());
  deallog << "cmpl_vec.vmult(compl_vec, real_vec) True" << std::endl;

  cmpl_matrix.vmult_add(cmpl_vec, real_vec);
  Assert((cmpl_vec(0).real() == 4.0) && (cmpl_vec(0).imag() == 4.0),
         ExcInternalError());
  deallog << "cmpl_vec.vmult_add(cmpl_vec, real_vec) True" << std::endl;

  cmpl_vec.reinit(1);

  cmpl_matrix.Tvmult(cmpl_vec, real_vec);
  Assert((cmpl_vec(0).real() == 2.0) && (cmpl_vec(0).imag() == 2.0),
         ExcInternalError());
  deallog << "cmpl_vec.Tvmult(cmpl_vec, real_vec) True" << std::endl;

  cmpl_matrix.Tvmult_add(cmpl_vec, real_vec);
  Assert((cmpl_vec(0).real() == 4.0) && (cmpl_vec(0).imag() == 4.0),
         ExcInternalError());
  deallog << "cmpl_vec.Tvmult_add(cmpl_vec, real_vec) True" << std::endl;
}

int
main()
{
  initlog();

  deallog << std::setprecision(0) << std::fixed;

  check_operations<double>();
}
