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


// test LAPACKFullMatrix::trace() by comparing to FullMatrix

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

  deallog << "trace difference: " << (F.trace() - M.trace()) << std::endl;
}


int
main()
{
  initlog();
  deallog.get_file_stream().precision(3);

  const std::vector<unsigned int> sizes = {{17, 391}};
  for (const auto &s : sizes)
    {
      deallog << "size=" << s << std::endl;
      test<double>(s);
    }
}
