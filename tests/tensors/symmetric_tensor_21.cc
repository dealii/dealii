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


// check SymmetricTensor<4,dim>::operator= (double)

#include <deal.II/base/symmetric_tensor.h>

#include <deal.II/lac/vector.h>

#include "../tests.h"

int
main()
{
  initlog();
  deallog << std::setprecision(3);

  const unsigned int      dim = 3;
  SymmetricTensor<4, dim> t;
  t[0][0][0][0] = t[1][0][1][0] = t[1][1][1][1] = t[2][2][2][2] =
    t[2][0][2][0]                               = 3;

  deallog << t.norm() << std::endl;
  t = 0;
  deallog << t.norm() << std::endl;

  AssertThrow(t.norm() == 0, ExcInternalError());

  deallog << "OK" << std::endl;
}
