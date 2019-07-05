// ---------------------------------------------------------------------
//
// Copyright (C) 2016 by the deal.II authors
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

// Test that Point operations on a CUDA device can be used.

#include <deal.II/base/point.h>

#include "../tests.h"

template <int dim, typename Number>
__global__ void
miscellaneous_kernel()
{
  Point<dim, Number> p_1;
  Point<dim, Number> p_2(Tensor<1, dim, Number>{});
  if (dim == 1)
    Point<dim, Number> p(1.);
  if (dim == 2)
    Point<dim, Number> p(1., 2.);
  if (dim == 3)
    Point<dim, Number> p(1., 2., 3.);

  auto p_3 = Point<dim, Number>::unit_vector(0);

  auto entry_1 = p_1(0);
  p_1(0)       = Number{1.};

  auto p_4 = p_1 + Tensor<1, dim, Number>{};
  auto p_5 = p_1 - Tensor<1, dim, Number>{};
  auto t_1 = p_1 - p_2;
  auto p_6 = -p_3;
  auto p_7 = p_4 / 2.;
  auto p_8 = p_2 * 5.;

  auto s_1 = p_1 * t_1;
  auto s_2 = p_2.square();
  auto s_3 = p_3.distance(p_5);
  auto s_4 = p_4.distance_square(p_1);
}

template <int dim, typename Number>
void
test_gpu()
{
  // Miscellaneous
  miscellaneous_kernel<dim, Number><<<1, 1>>>();
  // Check that the kernel was launched correctly
  AssertCuda(cudaGetLastError());
  // Check that there was no problem during the execution of the kernel
  AssertCuda(cudaDeviceSynchronize());

  deallog << "OK" << std::endl;
}

int
main()
{
  initlog();

  init_cuda();

  test_gpu<1, double>();
  test_gpu<2, double>();
  test_gpu<3, float>();
  test_gpu<1, float>();
  test_gpu<2, float>();
  test_gpu<3, float>();
}
