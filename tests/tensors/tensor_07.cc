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


// check scalar_product(Tensor<2,dim>,Tensor<2,dim>)

#include <deal.II/base/tensor.h>

#include <deal.II/lac/vector.h>

#include "../tests.h"

template <int dim>
void
test_constant()
{
  Tensor<2, dim> t;
  for (unsigned int i = 0; i < dim; ++i)
    for (unsigned int j = 0; j < dim; ++j)
      t[i][j] = 2.;
  deallog << "Constant dim " << dim << '\t' << scalar_product(t, t)
          << " compare " << 4 * dim * dim << std::endl;
}


template <int dim>
void
test_equal()
{
  Tensor<2, dim> t;
  unsigned int   sum = 0;
  for (unsigned int i = 0; i < dim; ++i)
    for (unsigned int j = 0; j < dim; ++j)
      {
        t[i][j] = i + dim * j;
        sum += (i + dim * j) * (i + dim * j);
      }

  deallog << "Equal    dim " << dim << '\t' << scalar_product(t, t)
          << " compare " << sum << std::endl;
}


template <int dim>
void
test_unequal()
{
  Tensor<2, dim> s;
  Tensor<2, dim> t;
  unsigned int   sum = 0;
  for (unsigned int i = 0; i < dim; ++i)
    for (unsigned int j = 0; j < dim; ++j)
      {
        s[i][j] = i + dim * j;
        t[i][j] = dim * i + j;
        sum += (i + dim * j) * (dim * i + j);
      }

  deallog << "Unequal  dim " << dim << '\t' << scalar_product(s, t)
          << " compare " << sum << std::endl;
}


int
main()
{
  initlog();

  test_constant<2>();
  test_constant<3>();
  test_constant<4>();

  test_equal<2>();
  test_equal<3>();
  test_equal<4>();

  test_unequal<2>();
  test_unequal<3>();
  test_unequal<4>();
}
