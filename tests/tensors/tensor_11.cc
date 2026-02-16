// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


// check Tensor::operator/=(int) and Tensor::operator/(int)

#include <deal.II/base/tensor.h>

#include "../tests.h"

template <typename Number>
void
test()
{
  Number a[3][3] = {{1, 2, 3}, {3, 4, 5}, {6, 7, 8}};

  const unsigned int     dim = 3;
  Tensor<2, dim, Number> t(a);

  deallog << "Before: " << t << std::endl;
  Tensor<2, dim, Number> u = t / 2;
  deallog << "operator/: " << u << std::endl;
  t /= 2;
  deallog << "operator/=: " << t << std::endl;
}

int
main()
{
  initlog();
  deallog.push("Test double");
  test<double>();
  deallog.pop();
  deallog.push("Test int");
  test<int>();
  deallog.pop();
}
