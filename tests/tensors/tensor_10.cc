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


// check schur_product for Tensor<rank, dim> with rank = 0, 1, 2, 3 and dim = 2

#include <deal.II/base/tensor.h>

#include <deal.II/lac/vector.h>

#include "../tests.h"

template <int rank, int dim>
void
test_same(const Tensor<rank, dim> &t, const Tensor<rank, dim> &compare)
{
  deallog << "Constant rank " << rank << " and dim " << dim << '\t'
          << schur_product(t, t) << " compare " << compare << std::endl;
}

template <int rank, int dim>
void
test_different(const Tensor<rank, dim> &t,
               const Tensor<rank, dim> &t2,
               const Tensor<rank, dim> &compare)
{
  deallog << "Constant rank " << rank << " and dim " << dim << '\t'
          << schur_product(t, t2) << " compare " << compare << std::endl;
}

int
main()
{
  initlog();

  {
    Tensor<0, 2> a, b, compare;
    a       = 2;
    b       = 2;
    compare = 4;
    test_same(a, compare);
    test_different(a, b, compare);
  }

  {
    Tensor<1, 2> a, b, compare;
    for (size_t d = 0; d < 2; ++d)
      {
        a[d]       = 2;
        b[d]       = 2;
        compare[d] = 4;
      }
    test_same(a, compare);
    test_different(a, b, compare);
  }

  {
    Tensor<2, 2> a, b, compare;
    for (size_t c = 0; c < 2; ++c)
      for (size_t d = 0; d < 2; ++d)
        {
          a[c][d]       = 2;
          b[c][d]       = 2;
          compare[c][d] = 4;
        }
    test_same(a, compare);
    test_different(a, b, compare);
  }

  {
    Tensor<3, 2> a, b, compare;

    for (size_t c = 0; c < 2; ++c)
      for (size_t d = 0; d < 2; ++d)
        for (size_t e = 0; e < 2; ++e)
          {
            a[c][d][e]       = 2;
            b[c][d][e]       = 2;
            compare[c][d][e] = 4;
          }
    test_same(a, compare);
    test_different(a, b, compare);
  }

  return 0;
}
