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


// check
// Tensor<1,dim,std::complex<double>>::operator*(Tensor<1,dim,std::complex<double>
// >) These tests are simple enough to do by hand.

#include <deal.II/base/tensor.h>

#include "../tests.h"

template <int dim>
void
test_tensor_01()
{
  Tensor<1, dim, std::complex<double>> t;
  for (unsigned int i = 0; i < dim; ++i)
    {
      t[i] = std::complex<double>(i, 1);
    }
  const Tensor<1, dim, std::complex<double>> res =
    t * std::complex<double>(1, 2);

  deallog << "dim = " << dim << ": " << res << std::endl;
}


template <int dim>
void
test_tensor_02()
{
  Tensor<1, dim, std::complex<double>> t;
  for (unsigned int i = 0; i < dim; ++i)
    {
      t[i] = std::complex<double>(1, i);
    }
  const Tensor<1, dim, std::complex<double>> res =
    t * std::complex<double>(1, 2);

  deallog << "dim = " << dim << ": " << res << std::endl;
}

int
main()
{
  initlog();
  deallog << std::setprecision(3);

  deallog << "First test" << std::endl;

  test_tensor_01<1>();
  test_tensor_01<2>();
  test_tensor_01<3>();

  deallog << "Second test" << std::endl;

  test_tensor_02<1>();
  test_tensor_02<2>();
  test_tensor_02<3>();

  deallog << "OK" << std::endl;
}
