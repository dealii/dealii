// ---------------------------------------------------------------------
//
// Copyright (C) 2017 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#include <array>
#include <deal.II/base/exceptions.h>
#include <iostream>
#include <numeric>

using namespace dealii;

__global__ void
double_value(double* x, double* y)
{
  y[threadIdx.x] = 2. * x[threadIdx.x];
}

int
main()
{
  constexpr int n = 4;

  std::array<double, n> host_x{};
  std::iota(host_x.begin(), host_x.end(), 1);
  std::array<double, n> host_y{};

  // Copy input data to device.
  double*     device_x;
  double*     device_y;
  cudaError_t cuda_error = cudaMalloc(&device_x, n * sizeof(double));
  AssertCuda(cuda_error);
  cuda_error = cudaMalloc(&device_y, n * sizeof(double));
  AssertCuda(cuda_error);
  cuda_error = cudaMemcpy(
    device_x, host_x.data(), n * sizeof(double), cudaMemcpyHostToDevice);
  AssertCuda(cuda_error);

  // Launch the kernel.
  double_value<<<1, n>>>(device_x, device_y);

  // Copy output data to host.
  cuda_error = cudaDeviceSynchronize();
  AssertCuda(cuda_error);
  cuda_error = cudaMemcpy(
    host_y.data(), device_y, n * sizeof(double), cudaMemcpyDeviceToHost);
  AssertCuda(cuda_error);

  // Print the results and test
  for(int i = 0; i < n; ++i)
    {
      std::cout << "y[" << i << "] = " << host_y[i] << "\n";
      AssertThrow(std::abs(host_y[i] - 2 * host_x[i]) < 1.e-10,
                  ExcInternalError());
    }

  cuda_error = cudaDeviceReset();
  AssertCuda(cuda_error);
  return 0;
}
