// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 1998 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Same as tensor.cc, but uses tensors based on std::complex<double> instead
// of double

#include <deal.II/base/tensor.h>

#include <deal.II/lac/vector.h>

#include <complex>

#include "../tests.h"

int
main()
{
  initlog();
  deallog << std::setprecision(3);

  double a_double[3][3][2] = {{{1, -1}, {2, 0}, {3, 0}},
                              {{3, 0}, {4, 0}, {5, 0}},
                              {{6, 0}, {7, 0}, {8, 3}}};
  double b_double[3][3][2] = {{{24, -2}, {31, -2}, {37, 6}},
                              {{45, -3}, {57, 0}, {69, 15}},
                              {{75, 12}, {96, 21}, {108, 48}}};

  const unsigned int   dim = 3;
  std::complex<double> a[dim][dim], b[dim][dim];
  for (unsigned int d = 0; d < dim; ++d)
    for (unsigned int e = 0; e < dim; ++e)
      {
        a[d][e] = std::complex<double>(a_double[d][e][0], a_double[d][e][1]);
        b[d][e] = std::complex<double>(b_double[d][e][0], b_double[d][e][1]);
      }

  Tensor<2, dim, std::complex<double>> t(a);
  Tensor<2, dim, std::complex<double>> tt;
  Tensor<2, dim, std::complex<double>> result(b);
  AssertThrow(transpose(transpose(t)) == t, ExcInternalError());
  AssertThrow(transpose(transpose(result)) == result, ExcInternalError());

  Vector<std::complex<double>> unrolled(9);

  // cast result to double to profit from zero
  // threshold and so on
  t.unroll(unrolled.data(), unrolled.data() + 9);
  deallog << "unrolled:";
  for (unsigned i = 0; i < 9; ++i)
    deallog << ' ' << unrolled(i);
  deallog << std::endl;

  deallog << "t=" << std::endl;
  for (unsigned int i = 0; i < dim; ++i)
    {
      for (unsigned int j = 0; j < dim; ++j)
        deallog << t[i][j] << ' ';
      deallog << std::endl;
    };

  deallog << "norm(t)=" << t.norm() << std::endl;

  tt = t * t;

  deallog << "tt=" << std::endl;
  for (unsigned int i = 0; i < dim; ++i)
    {
      for (unsigned int j = 0; j < dim; ++j)
        deallog << tt[i][j] << ' ';
      deallog << std::endl;
    };

  if (true)
    {
      deallog.push("Cross product");
      Tensor<1, 3, std::complex<double>> e1;
      Tensor<1, 3, std::complex<double>> e2;
      Tensor<1, 3, std::complex<double>> e3;
      e1[0]                                     = 1.;
      e2[1]                                     = 1.;
      e3[2]                                     = 1.;
      Tensor<1, 3, std::complex<double>> result = cross_product_3d(e1, e2);
      deallog << '\t' << result[0] << '\t' << result[1] << '\t' << result[2]
              << std::endl;

      result = cross_product_3d(e2, e3);
      deallog << '\t' << result[0] << '\t' << result[1] << '\t' << result[2]
              << std::endl;

      result = cross_product_3d(e3, e1);
      deallog << '\t' << result[0] << '\t' << result[1] << '\t' << result[2]
              << std::endl;

      deallog.pop();
    }

  if (tt == result)
    {
      deallog << "Result OK." << std::endl;
    }
  else
    {
      deallog << "Result WRONG!" << std::endl;
    };

  t = 0;
  deallog << t << std::endl;
}
