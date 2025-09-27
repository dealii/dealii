// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test LAPACKFullMatrix::compute_cholesky_factorization by comparing with
// FullMatrix

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/lapack_full_matrix.h>
#include <deal.II/lac/vector.h>

#include <iostream>

#include "../tests.h"

#include "create_matrix.h"


template <typename NumberType>
void
test(const unsigned int size)
{
  // Full matrix:
  FullMatrix<NumberType> F(size), C(size);
  create_spd(F);
  C.cholesky(F);

  // Lapack:
  LAPACKFullMatrix<NumberType> M(size);
  M = F;
  M.set_property(LAPACKSupport::symmetric);
  M.compute_cholesky_factorization();
  // factorization is stored in the lower diagonal part
  for (unsigned int i = 0; i < size; ++i)
    for (unsigned int j = i + 1; j < size; ++j)
      M(i, j) = 0.;

  FullMatrix<NumberType> diff(size);
  diff = M;
  diff.add(-1., C);

  const NumberType error = diff.frobenius_norm();
  deallog << size << " : " << diff.frobenius_norm() << std::endl;
  if (false)
    {
      std::cout << "Lapack:" << std::endl;
      M.print_formatted(std::cout);
      std::cout << "FullMatrix:" << std::endl;
      C.print_formatted(std::cout);
      std::cout << std::flush;
      AssertThrow(false, ExcInternalError());
    }
}


int
main()
{
  initlog();
  deallog.get_file_stream().precision(3);

  const std::vector<unsigned int> sizes = {{1, 3, 11, 17, 32, 64, 200, 391}};
  for (const auto &s : sizes)
    {
      // test<float>(s);
      test<double>(s);
    }
}
