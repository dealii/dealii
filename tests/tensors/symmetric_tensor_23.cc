// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2018 by the deal.II authors
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


// check operator<< for SymmetricTensor<2,dim> and SymmetricTensor<4,dim>

#include <deal.II/base/symmetric_tensor.h>

#include <deal.II/lac/vector.h>

#include "../tests.h"

int
main()
{
  initlog();
  deallog << std::setprecision(3);

  {
    SymmetricTensor<2, 1> t;
    t[0][0] = 1;

    double x[1] = {1};
    AssertThrow((t == SymmetricTensor<2, 1>(x)), ExcInternalError());
  }

  {
    SymmetricTensor<2, 2> t;
    t[0][0] = 1;
    t[1][1] = 2;
    t[0][1] = 3;

    double x[3] = {1, 2, 3};
    AssertThrow((t == SymmetricTensor<2, 2>(x)), ExcInternalError());
  }

  {
    SymmetricTensor<2, 3> t;
    t[0][0] = 1;
    t[1][1] = 2;
    t[2][2] = 3;
    t[0][1] = 4;
    t[0][2] = 5;
    t[1][2] = 6;

    double x[6] = {1, 2, 3, 4, 5, 6};
    AssertThrow((t == SymmetricTensor<2, 3>(x)), ExcInternalError());
  }

  deallog << "OK" << std::endl;
}
