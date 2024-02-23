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


// Test Givens rotations.


#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/utilities.h>
#include <deal.II/lac/vector.h>

#include "../tests.h"

template <typename NumberType>
void
test(const NumberType a, const NumberType b)
{
  FullMatrix<NumberType> rotation(2);
  Vector<NumberType>     x(2), y(2), res(2);

  x[0] = a;
  x[1] = b;
  y[1] = NumberType();

  const std::array<NumberType, 3> csr =
    Utilities::LinearAlgebra::givens_rotation(a, b);

  rotation(0, 0) = csr[0];  //  c
  rotation(1, 1) = csr[0];  //  c
  rotation(0, 1) = csr[1];  //  s
  rotation(1, 0) = -csr[1]; // -s
  y[0]           = csr[2];  //  r

  rotation.residual(res, x, y);

  const double norm = res.l2_norm();
  deallog << norm << std::endl;

  if (norm > 1e-12 || csr[2] < 0.)
    {
      deallog << "x:" << std::endl;
      x.print(deallog.get_file_stream());
      deallog << "Givens:" << std::endl;
      rotation.print(deallog.get_file_stream(), 10, 6);
      deallog << "y:" << std::endl;
      y.print(deallog.get_file_stream());
      deallog << "res:" << std::endl;
      res.print(deallog.get_file_stream());
      AssertThrow(false, ExcInternalError());
    }
}

int
main()
{
  initlog();
  deallog << std::setprecision(6);

  deallog << "Residuals:" << std::endl;
  // g==0
  test<double>(3., 0.);
  test<double>(-2., 0.);
  // f==0
  test<double>(0., 2.);
  test<double>(0., -5.);
  // |f| > |g|
  test<double>(15., 3.);
  test<double>(15., -4.);
  test<double>(-17., 2.);
  test<double>(-18., -5.);
  // |f| < |g|
  test<double>(2., -4.);
  test<double>(-2., 3.);
  test<double>(-3., -7.);
  test<double>(5., 9.);
}
