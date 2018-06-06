// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2017 by the deal.II authors
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


// test full 2x2 tensors

#include <deal.II/base/tensor.h>

#include "../tests.h"

int
main()
{
  std::ofstream logfile("output");
  deallog << std::setprecision(3);
  deallog.attach(logfile);

  Tensor<2, 2> t;
  t[0][0] = 1;
  t[1][1] = 2;
  t[0][1] = 4;
  t[1][0] = 4;
  // make sure transposition doesn't change
  // anything
  AssertThrow(t == transpose(t), ExcInternalError());

  // check norm of tensor
  AssertThrow(std::fabs(t.norm() - std::sqrt(1. * 1 + 2 * 2 + 2 * 4 * 4)) <
                1e-14,
              ExcInternalError());

  deallog << "OK" << std::endl;
}
