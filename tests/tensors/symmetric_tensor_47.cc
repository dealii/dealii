// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// check  SymmetricTensor::symmetrize on fourth order Tensors


#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/vectorization.h>

#include "../tests.h"


template <int dim, typename Number>
void
make_fourth_order_tensor(Tensor<4, dim, Number> &t)
{
  for (unsigned int i = 0; i < dim; ++i)
    for (unsigned int j = 0; j < dim; ++j)
      for (unsigned int k = 0; k < dim; ++k)
        for (unsigned int l = 0; l < dim; ++l)
          t[i][j][k][l] = (i + 1) * (j + 2) * (k + 3) * (l + 4);
}

template <int dim, typename Number>
void
print_fourth_order_tensor(const SymmetricTensor<4, dim, Number> &t)
{
  const Number tolerance = 100. * std::numeric_limits<Number>::epsilon();
  for (unsigned int i = 0; i < dim; ++i)
    for (unsigned int j = 0; j < dim; ++j)
      for (unsigned int k = 0; k < dim; ++k)
        for (unsigned int l = 0; l < dim; ++l)
          deallog << i << ' ' << j << ' ' << k << ' ' << l << ' '
                  << filter_out_small_numbers(t[i][j][k][l], tolerance)
                  << std::endl;
}


int
main()
{
  initlog();
  const unsigned int     dim = 3;
  Tensor<4, dim, double> t;
  make_fourth_order_tensor(t);

  SymmetricTensor<4, dim, double> st_minor;

  st_minor = symmetrize(t, false);

  // test minor symmetry
  for (unsigned int i = 0; i < dim; ++i)
    for (unsigned int j = 0; j < dim; ++j)
      for (unsigned int k = 0; k < dim; ++k)
        for (unsigned int l = 0; l < dim; ++l)
          AssertThrow(st_minor[i][j][k][l] == st_minor[i][j][l][k] &&
                        st_minor[j][i][k][l] == st_minor[i][j][k][l] &&
                        st_minor[i][j][k][l] == st_minor[j][i][l][k],
                      ExcInternalError());

  print_fourth_order_tensor(st_minor);

  deallog << "MINOR SYMMETRY OK" << std::endl;


  SymmetricTensor<4, dim, double> st_major;

  st_major = symmetrize(t, true);

  // test major symmetry
  for (unsigned int i = 0; i < dim; ++i)
    for (unsigned int j = 0; j < dim; ++j)
      for (unsigned int k = 0; k < dim; ++k)
        for (unsigned int l = 0; l < dim; ++l)
          AssertThrow(st_major[i][j][k][l] == st_major[k][l][i][j],
                      ExcInternalError());

  print_fourth_order_tensor(st_major);

  deallog << "MAJOR SYMMETRY OK" << std::endl;
  deallog << "OK" << std::endl;
}
