// ---------------------------------------------------------------------
//
// Copyright (C) 2014 - 2020 by the deal.II authors
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


// Tests LAPACKFullMatrix::Tmmult with  A==B which shall use Xsyrk instead Xgemm

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/lapack_full_matrix.h>
#include <deal.II/lac/vector.h>

#include <iostream>
#include <tuple>

#include "../tests.h"

#include "create_matrix.h"

DeclException5(ExcEl,
               int,
               int,
               double,
               double,
               double,
               << "Error in element (" << arg1 << "," << arg2 << "): " << arg3
               << " != " << arg4 << " delta=" << arg5);

template <typename NumberType>
void
test(const unsigned int n, const unsigned int k, const NumberType eps)
{
  deallog << n << " " << k << " " << std::endl;
  FullMatrix<NumberType>       A(n, k), C(n, n);
  LAPACKFullMatrix<NumberType> AL(n, k), CL(n, n);

  create_random(AL);

  A = AL;

  A.mTmult(C, A);
  AL.mTmult(CL, AL);

  for (unsigned int i = 0; i < n; ++i)
    for (unsigned int j = 0; j < n; ++j)
      AssertThrow(std::abs(C(i, j) - CL(i, j)) < eps * std::abs(CL(i, j)),
                  ExcEl(i, j, C(i, j), CL(i, j), C(i, j) - CL(i, j)));

  deallog << "OK" << std::endl;

  AL.mTmult(CL, AL, true);

  for (unsigned int i = 0; i < n; ++i)
    for (unsigned int j = 0; j < n; ++j)
      AssertThrow(std::abs(2. * C(i, j) - CL(i, j)) < eps * std::abs(CL(i, j)),
                  ExcEl(i, j, 2. * C(i, j), CL(i, j), 2. * C(i, j) - CL(i, j)));

  deallog << "OK adding" << std::endl;
}

int
main()
{
  initlog();
  deallog.get_file_stream().precision(3);

  const std::vector<std::array<unsigned int, 2>> sizes = {{{3, 3}},
                                                          {{7, 7}},
                                                          {{51, 51}},
                                                          {{320, 320}},
                                                          {{3, 9}},
                                                          {{9, 7}},
                                                          {{10, 5}},
                                                          {{320, 120}}};

  deallog.push("double");
  for (auto el : sizes)
    test<double>(el[0], el[1], 1e-13);
  deallog.pop();

  deallog.push("float");
  for (auto el : sizes)
    test<float>(el[0], el[1], 1e-5);
  deallog.pop();
}
