// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test for mixed Number type operations of Tensor<0, dim, Number> and
// Tensor<1, dim. Number> (addition, subtraction, contraction).

#include <deal.II/base/tensor.h>

#include <complex>

#include "../tests.h"

int
main()
{
  initlog();

  {
    dealii::Tensor<0, 3, float>                f = 2.;
    dealii::Tensor<0, 3, double>               d = 4.;
    dealii::Tensor<0, 3, std::complex<double>> c = 8.;

    deallog << f + d << std::endl;
    deallog << f - d << std::endl;
    deallog << f * d << std::endl;

    deallog << d + c << std::endl;
    deallog << d - c << std::endl;
    deallog << d * c << std::endl;

    float f_scalar = 10.;
    deallog << (d * f_scalar) << std::endl;
    deallog << (d / f_scalar) << std::endl;
    deallog << (f_scalar * d) << std::endl;
    deallog << (d *= f_scalar) << std::endl;
    deallog << (d /= f_scalar) << std::endl;

    double d_scalar = 10.;
    deallog << (c * d_scalar) << std::endl;
    deallog << (c / d_scalar) << std::endl;
    deallog << (d_scalar * c) << std::endl;
    deallog << (c *= d_scalar) << std::endl;
    deallog << (c /= d_scalar) << std::endl;
  }

  {
    dealii::Tensor<1, 2, float>                f({1., 2.});
    dealii::Tensor<1, 2, double>               d({4., 8.});
    dealii::Tensor<1, 2, std::complex<double>> c({16., 32.});

    deallog << f + d << std::endl;
    deallog << f - d << std::endl;
    deallog << f * d << std::endl;

    deallog << d + c << std::endl;
    deallog << d - c << std::endl;
    deallog << d * c << std::endl;

    float f_scalar = 10.;
    deallog << (d * f_scalar) << std::endl;
    deallog << (d / f_scalar) << std::endl;
    deallog << (f_scalar * d) << std::endl;
    deallog << (d *= f_scalar) << std::endl;
    deallog << (d /= f_scalar) << std::endl;

    double d_scalar = 10.;
    deallog << (c * d_scalar) << std::endl;
    deallog << (c / d_scalar) << std::endl;
    deallog << (d_scalar * c) << std::endl;
    deallog << (c *= d_scalar) << std::endl;
    deallog << (c /= d_scalar) << std::endl;
  }

  {
    dealii::Tensor<1, 3, float>                f({1., 2., 4.});
    dealii::Tensor<1, 3, double>               d({4., 8., 16.});
    dealii::Tensor<1, 3, std::complex<double>> c({32., 64., 128.});

    deallog << f + d << std::endl;
    deallog << f - d << std::endl;
    deallog << f * d << std::endl;

    deallog << d + c << std::endl;
    deallog << d - c << std::endl;
    deallog << d * c << std::endl;

    float f_scalar = 10.;
    deallog << (d * f_scalar) << std::endl;
    deallog << (d / f_scalar) << std::endl;
    deallog << (f_scalar * d) << std::endl;
    deallog << (d *= f_scalar) << std::endl;
    deallog << (d /= f_scalar) << std::endl;

    double d_scalar = 10.;
    deallog << (c * d_scalar) << std::endl;
    deallog << (c / d_scalar) << std::endl;
    deallog << (d_scalar * c) << std::endl;
    deallog << (c *= d_scalar) << std::endl;
    deallog << (c /= d_scalar) << std::endl;
  }

  deallog << "OK!" << std::endl;
}
