// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2018 by the deal.II authors
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


// test that the product between a tensor and an integer behaves just
// like that of the tensor and the integer-converted-to-double

#include <deal.II/base/template_constraints.h>
#include <deal.II/base/tensor.h>

#include <complex>
#include <typeinfo>

#include "../tests.h"



int
main()
{
  initlog();

  // try it for a rank-1 tensor
  {
    Tensor<1, 3> t;
    t[0] = 1.23456;
    t[1] = 2.46802;
    t[2] = 3.69258;
    AssertThrow(7 * t == 7.0 * t, ExcInternalError());
    AssertThrow(t * 7 == t * 7.0, ExcInternalError());
    AssertThrow(t * 7 == 7 * t, ExcInternalError());
    AssertThrow((t * 7 - (t + t + t + t + t + t + t)).norm() < 1e-12,
                ExcInternalError());
  }

  // now also try it for a rank-2 tensor (higher rank tensors are
  // composed of lower-rank tensors)
  {
    Tensor<2, 2> t;
    t[0][0] = 1.23456;
    t[0][1] = 7.87965;
    t[1][0] = 2.64686;
    t[1][1] = 3.35792;
    AssertThrow(7 * t == 7.0 * t, ExcInternalError());
    AssertThrow(t * 7 == t * 7.0, ExcInternalError());
    AssertThrow(t * 7 == 7 * t, ExcInternalError());
    AssertThrow((t * 7 - (t + t + t + t + t + t + t)).norm() < 1e-12,
                ExcInternalError());
  }

  deallog << "OK" << std::endl;
}
