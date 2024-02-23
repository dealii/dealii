// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2014 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test LAPACKFullMatrix:: norms for non-symmetric matrices

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
  FullMatrix<NumberType> F(size);
  create_random(F);

  // Lapack:
  LAPACKFullMatrix<NumberType> M(size);
  M = F;

  deallog << "l1:        " << (F.l1_norm() - M.l1_norm()) << std::endl
          << "linfty:    " << (F.linfty_norm() - M.linfty_norm()) << std::endl
          << "frobenius: " << (F.frobenius_norm() - M.frobenius_norm())
          << std::endl;
}


int
main()
{
  initlog();
  deallog.get_file_stream().precision(3);

  const std::vector<unsigned int> sizes = {{1, 3, 11, 17, 32, 64, 200, 391}};
  for (const auto &s : sizes)
    {
      deallog << "size=" << s << std::endl;
      // test<float>(s);
      test<double>(s);
    }
}
