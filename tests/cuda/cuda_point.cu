// ---------------------------------------------------------------------
//
// Copyright (C) 2019 by the deal.II authors
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
miscellaneous_kernel(Number check[16])
{
  Point<dim, Number> p_1;
  check[0] = p_1.norm_square();
  Point<dim, Number> p_2(Tensor<1, dim, Number>{});
  check[1] = p_2.norm_square();
  if (dim == 1)
    {
      Point<dim, Number> p(1.);
      check[2] = p.norm_square();
    }
  if (dim == 2)
    {
      Point<dim, Number> p(.6, .8);
      check[2] = p.norm_square();
    }
  if (dim == 3)
    {
      Point<dim, Number> p(.48, .64, .6);
      check[2] = p.norm_square();
    }

  auto p_3 = Point<dim, Number>::unit_vector(0);
  check[3] = p_3.norm_square();

  auto entry_1 = p_1(0);
  check[4]     = entry_1;
  p_1(0)       = Number{1.};
  check[5]     = p_1.norm_square();
  auto p_4     = p_1 + Tensor<1, dim, Number>{};
  check[6]     = p_4.norm_square();
  auto p_5     = p_1 - Tensor<1, dim, Number>{};
  check[7]     = p_5.norm_square();
  auto t_1     = p_1 - p_2;
  check[8]     = t_1.norm_square();
  auto p_6     = -p_3;
  check[9]     = p_6.norm_square();
  auto p_7     = p_4 / 2.;
  check[10]    = p_7.norm_square();
  auto p_8     = p_7 * 5.;
  check[11]    = p_8.norm_square();

  auto s_1  = p_1 * t_1;
  check[12] = s_1;
  auto s_2  = p_2.square();
  check[13] = s_2;
  auto s_3  = p_3.distance(p_5);
  check[14] = s_3;
  auto s_4  = p_4.distance_square(p_1);
  check[15] = s_4;
}

template <int dim, typename Number>
void
test_gpu()
{
  Number *           check;
  const unsigned int n_tests = 16;

  auto cuda_error = cudaMalloc(&check, n_tests * sizeof(Number));
  AssertCuda(cuda_error);

  // Miscellaneous
  miscellaneous_kernel<dim, Number><<<1, 1>>>(check);
  // Check that the kernel was launched correctly
  AssertCuda(cudaPeekAtLastError());
  // Check that there was no problem during the execution of the kernel
  AssertCuda(cudaDeviceSynchronize());

  std::vector<Number> check_host(n_tests);

  cuda_error = cudaMemcpy(check_host.data(),
                          check,
                          n_tests * sizeof(Number),
                          cudaMemcpyDeviceToHost);
  AssertCuda(cuda_error);

  const double tolerance = 1.e-8;
  AssertThrow(std::abs(check_host[0] - 0.) < tolerance, ExcInternalError());
  AssertThrow(std::abs(check_host[1] - 0.) < tolerance, ExcInternalError());
  AssertThrow(std::abs(check_host[2] - 1.) < tolerance, ExcInternalError());
  AssertThrow(std::abs(check_host[3] - 1.) < tolerance, ExcInternalError());
  AssertThrow(std::abs(check_host[4] - 0.) < tolerance, ExcInternalError());
  AssertThrow(std::abs(check_host[5] - 1.) < tolerance, ExcInternalError());
  AssertThrow(std::abs(check_host[6] - 1.) < tolerance, ExcInternalError());
  AssertThrow(std::abs(check_host[7] - 1.) < tolerance, ExcInternalError());
  AssertThrow(std::abs(check_host[8] - 1.) < tolerance, ExcInternalError());
  AssertThrow(std::abs(check_host[9] - 1.) < tolerance, ExcInternalError());
  AssertThrow(std::abs(check_host[10] - .25) < tolerance, ExcInternalError());
  AssertThrow(std::abs(check_host[11] - 6.25) < tolerance, ExcInternalError());
  AssertThrow(std::abs(check_host[12] - 1.) < tolerance, ExcInternalError());
  AssertThrow(std::abs(check_host[13] - 0.) < tolerance, ExcInternalError());
  AssertThrow(std::abs(check_host[14] - 0.) < tolerance, ExcInternalError());
  AssertThrow(std::abs(check_host[15] - 0.) < tolerance, ExcInternalError());

  cuda_error = cudaFree(check);
  AssertCuda(cuda_error);

  deallog << "OK" << std::endl;
}

int
main()
{
  initlog();

  init_cuda();

  test_gpu<1, double>();
  test_gpu<2, double>();
  test_gpu<3, double>();
  test_gpu<1, float>();
  test_gpu<2, float>();
  test_gpu<3, float>();
}
