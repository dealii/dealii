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


// Same as tensor.cc, but uses tensors based on floats instead of doubles

#include <deal.II/base/tensor.h>

#include <deal.II/lac/vector.h>

#include "../tests.h"

int
main()
{
  initlog();
  deallog << std::setprecision(3);

  float a[3][3] = {{1, 2, 3}, {3, 4, 5}, {6, 7, 8}};
  float b[3][3] = {{25, 31, 37}, {45, 57, 69}, {75, 96, 117}};

  const unsigned int    dim = 3;
  Tensor<2, dim, float> t(a);
  Tensor<2, dim, float> tt;
  Tensor<2, dim, float> result(b);
  AssertThrow(transpose(transpose(t)) == t, ExcInternalError());
  AssertThrow(transpose(transpose(result)) == result, ExcInternalError());

  Vector<float> unrolled(9);

  // cast result to double to profit from zero
  // threshold and so on
  t.unroll(unrolled.data(), unrolled.data() + 9);
  deallog << "unrolled:";
  for (unsigned i = 0; i < 9; ++i)
    deallog << ' ' << static_cast<double>(unrolled(i));
  deallog << std::endl;

  deallog << "t=" << std::endl;
  for (unsigned int i = 0; i < dim; ++i)
    {
      for (unsigned int j = 0; j < dim; ++j)
        deallog << static_cast<double>(t[i][j]) << ' ';
      deallog << std::endl;
    };

  deallog << "norm(t)=" << t.norm() << std::endl;

  tt = t * t;

  deallog << "tt=" << std::endl;
  for (unsigned int i = 0; i < dim; ++i)
    {
      for (unsigned int j = 0; j < dim; ++j)
        deallog << static_cast<double>(tt[i][j]) << ' ';
      deallog << std::endl;
    };

  if (true)
    {
      deallog.push("Cross product");
      Tensor<1, 3, float> e1;
      Tensor<1, 3, float> e2;
      Tensor<1, 3, float> e3;
      e1[0]                      = 1.;
      e2[1]                      = 1.;
      e3[2]                      = 1.;
      Tensor<1, 3, float> result = cross_product_3d(e1, e2);
      deallog << '\t' << static_cast<double>(result[0]) << '\t'
              << static_cast<double>(result[1]) << '\t'
              << static_cast<double>(result[2]) << std::endl;

      result = cross_product_3d(e2, e3);
      deallog << '\t' << static_cast<double>(result[0]) << '\t'
              << static_cast<double>(result[1]) << '\t'
              << static_cast<double>(result[2]) << std::endl;

      result = cross_product_3d(e3, e1);
      deallog << '\t' << static_cast<double>(result[0]) << '\t'
              << static_cast<double>(result[1]) << '\t'
              << static_cast<double>(result[2]) << std::endl;

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
}
