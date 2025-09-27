// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Based on tensor_complex.cc. It tests cross_product_3d when one tensor is
// std::complex<double> and the other is double.

#include <deal.II/base/tensor.h>

#include <deal.II/lac/vector.h>

#include <complex>

#include "../tests.h"

int
main()
{
  initlog();
  deallog << std::setprecision(3);

  deallog.push("Cross product");
  Tensor<1, 3, std::complex<double>> e1;
  Tensor<1, 3, double>               e2;
  Tensor<1, 3, std::complex<double>> e3;
  e1[0]                                     = std::complex<double>(1., 1.);
  e2[1]                                     = 1.;
  e3[2]                                     = std::complex<double>(1., -1.);
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
