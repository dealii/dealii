// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2020 by the deal.II authors
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

// Test that we can convert a SymmetricTensor<4,dim> to a Tensor<4,dim>

#include <deal.II/base/symmetric_tensor.h>

#include "../tests.h"


template <int dim>
void initialize(SymmetricTensor<2, dim> &st)
{
  for (unsigned int i = 0; i < dim; ++i)
    for (unsigned int j = i; j < dim; ++j)
      st[i][j] = (i + 1) * dim + (j - i);
}


template <int dim>
void initialize(SymmetricTensor<4, dim> &st)
{
  for (unsigned int i = 0; i < dim; ++i)
    for (unsigned int j = i; j < dim; ++j)
      for (unsigned int k = 0; k < dim; ++k)
        for (unsigned int l = k; l < dim; ++l)
          st[i][j][k][l] = (i + 1) * dim + (j - i) + (k + 1) * dim + (l - k);
}



template <int rank, int dim>
void
check()
{
  // build a regular tensor
  SymmetricTensor<rank, dim> st;
  initialize(st);
  deallog << "st=" << st << std::endl;

  Tensor<rank, dim> t(st); // check conversion constructor
  deallog << "t =" << t << std::endl;

  t = st; // check assignment operator
  deallog << "t =" << t << std::endl;
}


int
main()
{
  initlog();
  deallog << std::setprecision(3);

  deallog << "checking rank 2 tensors" << std::endl;
  check<2, 1>();
  check<2, 2>();
  check<2, 3>();

  deallog << "checking rank 4 tensors" << std::endl;
  check<4, 1>();
  check<4, 2>();
  check<4, 3>();

  deallog << "OK" << std::endl;
}
