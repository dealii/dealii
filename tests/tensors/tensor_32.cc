// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2020 by the deal.II authors
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


// Test that we can initialize a Tensor from ArrayView

#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/tensor.h>

#include "../tests.h"


template <int dim>
void
test()
{
  const std::vector<double> values = {1, 2, 3, 4, 5, 6, 7, 8, 9};

  // Check for Tensor<1,dim>
  deallog << Tensor<1, dim>(
               make_array_view(values.begin(), values.begin() + dim))
          << std::endl;

  // Check for Tensor<2,dim>
  deallog << Tensor<2, dim>(
               make_array_view(values.begin(), values.begin() + dim * dim))
          << std::endl;
}

int
main(int argc, char *argv[])
{
  initlog();

  test<1>();
  test<2>();
  test<3>();

  deallog << "OK" << std::endl;

  return 0;
}
