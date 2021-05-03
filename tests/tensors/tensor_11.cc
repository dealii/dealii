// ---------------------------------------------------------------------
//
// Copyright (C) 2019 by the deal.II authors
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
