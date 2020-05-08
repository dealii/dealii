// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2020 by the deal.II authors
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
