// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2020 by the deal.II authors
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


// check Tensor::operator= (double)

#include <deal.II/base/tensor.h>

#include <deal.II/lac/vector.h>

#include "../tests.h"

int
main()
{
  initlog();
  deallog << std::setprecision(3);

  double a[3][3] = {{1, 2, 3}, {3, 4, 5}, {6, 7, 8}};

  const unsigned int dim = 3;
  Tensor<2, dim>     t(a);

  deallog << t.norm() << std::endl;
  t = 0;
  deallog << t.norm() << std::endl;

  Assert(t.norm() == 0, ExcInternalError());

  deallog << "OK" << std::endl;
}
