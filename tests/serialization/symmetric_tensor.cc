// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2010 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


// check serialization for SymmetricTensor<2,dim>

#include <deal.II/base/symmetric_tensor.h>

#include "serialization.h"


void
test()
{
  const unsigned int dim  = 3;
  const unsigned int rank = 2;

  double a1[3][3] = {{1., 2., 3.}, {2., 5., 6.}, {3., 6., 9.}};
  SymmetricTensor<rank, dim> t1((Tensor<rank, dim>(a1)));


  double a2[3][3] = {{10., 11., 12.}, {11., 14., 15.}, {12., 15., 18.}};
  SymmetricTensor<rank, dim> t2((Tensor<rank, dim>(a2)));

  verify(t1, t2);
}


int
main()
{
  initlog();

  test();

  deallog << "OK" << std::endl;
}
