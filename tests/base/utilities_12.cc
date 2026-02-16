// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2017 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


// test Utilities::fixed_power on VectorizedArray, similar to utilities_02

#include <deal.II/base/utilities.h>
#include <deal.II/base/vectorization.h>

#include "../tests.h"


template <int dim>
void
test()
{
  VectorizedArray<double> v1 = make_vectorized_array(2.);
  v1                         = Utilities::fixed_power<dim>(v1);
  deallog << v1[0] << std::endl;
  v1 = -2;
  v1 = Utilities::fixed_power<dim>(v1);
  deallog << v1[0] << std::endl;
  v1 = 2.5;
  v1 = Utilities::fixed_power<dim>(v1);
  deallog << v1[0] << std::endl;
  VectorizedArray<float> v2 = make_vectorized_array<float>(-2.5);
  v2                        = Utilities::fixed_power<dim>(v2);
  deallog << (double)v2[0] << std::endl;
  deallog << std::endl;
}



int
main()
{
  initlog();

  test<0>();
  test<1>();
  test<2>();
  test<3>();
  test<4>();
}
