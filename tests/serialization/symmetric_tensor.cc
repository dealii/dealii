// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2018 by the deal.II authors
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
