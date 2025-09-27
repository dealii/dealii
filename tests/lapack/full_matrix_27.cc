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


// test LAPACKFullMatrix::rank1_update() for rank1 update of a matrix

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
  // Lapack:
  LAPACKFullMatrix<NumberType> A(size);
  A.set_property(LAPACKSupport::symmetric);
  Vector<NumberType> v(size);
  for (unsigned int i = 0; i < size; ++i)
    {
      v(i) = random_value<NumberType>();
      for (unsigned int j = i; j < size; ++j)
        {
          const NumberType val = random_value<NumberType>();
          A(i, j)              = val;
          A(j, i)              = val;
        }
    }


  const NumberType a = random_value<NumberType>();

  LAPACKFullMatrix<NumberType> B(size);
  B = A;

  B.rank1_update(a, v);

  for (unsigned int i = 0; i < size; ++i)
    for (unsigned int j = 0; j < size; ++j)
      {
        const NumberType diff = A(i, j) + a * v(i) * v(j) - B(i, j);
        AssertThrow(std::abs(diff) < 1e-10 * std::abs(B(i, j)),
                    ExcMessage("diff=" + std::to_string(diff)));
      }

  deallog << "OK" << std::endl;
}


int
main()
{
  initlog();
  deallog.get_file_stream().precision(3);

  const std::vector<unsigned int> sizes = {{17, 35, 391}};
  for (const auto &s : sizes)
    {
      deallog << "size=" << s << std::endl;
      test<double>(s);
    }
}
