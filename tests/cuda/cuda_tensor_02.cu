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

// Test operator[] and norm_square of cuda_tensor.

#include <deal.II/base/tensor.h>

#include "../tests.h"

template <int rank, int dim, typename Number>
__global__ void
miscellaneous_kernel()
{
  // constructors
  typename Tensor<rank, dim, Number>::array_type array{};
  Tensor<rank, dim, Number>                      dummy_1(array);
  Tensor<rank, dim, Number>                      dummy_2;
  Tensor<rank, dim, Number>                      dummy_3 = dummy_2;

  // access
  Tensor<rank + 1, dim, Number> initializer_1;
  const auto                    dummy_5 = initializer_1[0];

  // assignment
  dummy_2 = dummy_3;
}

template <int rank, int dim, typename Number>
__global__ void
summation_kernel(Tensor<rank, dim, Number> *t,
                 Tensor<rank, dim, Number> *t1,
                 Tensor<rank, dim, Number> *t2)
{
  *t2 += *t;
  *t1 = *t1 + *t;
}

template <int rank, int dim, typename Number>
__global__ void
subtraction_kernel(Tensor<rank, dim, Number> *t,
                   Tensor<rank, dim, Number> *t1,
                   Tensor<rank, dim, Number> *t2)
{
  *t2 -= *t;
  *t1 = *t1 - *t;
}

template <int rank, int dim, typename Number>
__global__ void
multiplication_kernel(Tensor<rank, dim, Number> *t,
                      Tensor<rank, dim, Number> *t1,
                      Tensor<rank, dim, Number> *t2)
{
  *t1 = *t * Number(2.);
  *t2 = Number(2.) * *t;
  *t *= 2.;
}

template <int rank, int dim, typename Number>
__global__ void
division_kernel(Tensor<rank, dim, Number> *t,
                Tensor<rank, dim, Number> *t1,
                Tensor<rank, dim, Number> *t2)
{
  *t1 = *t / Number(2.);
  *t /= 2.;
  *t2 = *t1;
}

template <int dim, typename Number>
__global__ void init_kernel(Tensor<0, dim, Number> *t)
{
  if (threadIdx.x == 0)
    *t = 1.;
}

template <int dim, typename Number>
__global__ void init_kernel(Tensor<1, dim, Number> *t)
{
  const unsigned int i = threadIdx.x;
  if (i < dim)
    (*t)[i] = i + 1.;
}

template <int dim, typename Number>
__global__ void init_kernel(Tensor<2, dim, Number> *t)
{
  const unsigned int i = threadIdx.y;
  const unsigned int j = threadIdx.x;
  if ((i < dim) && (j < dim))
    (*t)[i][j] = j + i * dim + 1.;
}


template <int rank, int dim, typename Number>
void
test_gpu()
{
  const double tolerance = 1.e-8;

  Tensor<rank, dim, Number> *t_dev;
  Tensor<rank, dim, Number> *t1_dev;
  Tensor<rank, dim, Number> *t2_dev;

  Tensor<rank, dim, Number> t_host;
  Tensor<rank, dim, Number> t1_host;
  Tensor<rank, dim, Number> t2_host;

  Tensor<rank, dim, Number> reference_host;

  // Allocate objects on the device
  cudaError_t cuda_error =
    cudaMalloc(&t_dev, sizeof(Tensor<rank, dim, Number>));
  AssertCuda(cuda_error);
  cuda_error = cudaMalloc(&t1_dev, sizeof(Tensor<rank, dim, Number>));
  AssertCuda(cuda_error);
  cuda_error = cudaMalloc(&t2_dev, sizeof(Tensor<rank, dim, Number>));
  AssertCuda(cuda_error);

  // Initialize
  dim3 block_dim(dim, dim);
  init_kernel<<<1, block_dim>>>(t_dev);
  cuda_error = cudaMemcpy(&reference_host,
                          t_dev,
                          sizeof(Tensor<rank, dim, Number>),
                          cudaMemcpyDeviceToHost);
  AssertCuda(cuda_error);

  // Test multiplication.
  multiplication_kernel<<<1, 1>>>(t_dev, t1_dev, t2_dev);

  cuda_error = cudaMemcpy(&t_host,
                          t_dev,
                          sizeof(Tensor<rank, dim, Number>),
                          cudaMemcpyDeviceToHost);
  AssertCuda(cuda_error);
  cuda_error = cudaMemcpy(&t1_host,
                          t1_dev,
                          sizeof(Tensor<rank, dim, Number>),
                          cudaMemcpyDeviceToHost);
  AssertCuda(cuda_error);
  cuda_error = cudaMemcpy(&t2_host,
                          t2_dev,
                          sizeof(Tensor<rank, dim, Number>),
                          cudaMemcpyDeviceToHost);
  AssertCuda(cuda_error);

  reference_host *= 2;
  AssertThrow((t_host - reference_host).norm() < tolerance, ExcInternalError());
  AssertThrow((t1_host - reference_host).norm() < tolerance,
              ExcInternalError());
  AssertThrow((t2_host - reference_host).norm() < tolerance,
              ExcInternalError());

  deallog << "multiplication OK" << std::endl;

  // Test division.
  division_kernel<<<1, 1>>>(t_dev, t1_dev, t2_dev);
  cuda_error = cudaMemcpy(&t_host,
                          t_dev,
                          sizeof(Tensor<rank, dim, Number>),
                          cudaMemcpyDeviceToHost);
  AssertCuda(cuda_error);
  cuda_error = cudaMemcpy(&t1_host,
                          t1_dev,
                          sizeof(Tensor<rank, dim, Number>),
                          cudaMemcpyDeviceToHost);
  AssertCuda(cuda_error);

  reference_host /= 2.;
  AssertThrow((t_host - reference_host).norm() < tolerance, ExcInternalError());
  AssertThrow((t1_host - reference_host).norm() < tolerance,
              ExcInternalError());

  deallog << "division OK" << std::endl;

  // Test summation
  summation_kernel<<<1, 1>>>(t_dev, t1_dev, t2_dev);
  cuda_error = cudaMemcpy(&t1_host,
                          t1_dev,
                          sizeof(Tensor<rank, dim, Number>),
                          cudaMemcpyDeviceToHost);
  AssertCuda(cuda_error);
  cuda_error = cudaMemcpy(&t2_host,
                          t2_dev,
                          sizeof(Tensor<rank, dim, Number>),
                          cudaMemcpyDeviceToHost);
  AssertCuda(cuda_error);

  reference_host *= 2.;
  AssertThrow((t1_host - reference_host).norm() < tolerance,
              ExcInternalError());
  AssertThrow((t2_host - reference_host).norm() < tolerance,
              ExcInternalError());


  // Test subtraction
  subtraction_kernel<<<1, 1>>>(t_dev, t1_dev, t2_dev);
  cuda_error = cudaMemcpy(&t1_host,
                          t1_dev,
                          sizeof(Tensor<rank, dim, Number>),
                          cudaMemcpyDeviceToHost);
  AssertCuda(cuda_error);
  cuda_error = cudaMemcpy(&t2_host,
                          t2_dev,
                          sizeof(Tensor<rank, dim, Number>),
                          cudaMemcpyDeviceToHost);

  reference_host /= 2.;
  AssertThrow((t1_host - reference_host).norm() < tolerance,
              ExcInternalError());
  AssertThrow((t2_host - reference_host).norm() < tolerance,
              ExcInternalError());

  // Miscellaneous
  miscellaneous_kernel<rank, dim, Number><<<1, 1>>>();

  // Free memory
  cuda_error = cudaFree(t_dev);
  AssertCuda(cuda_error);
  cuda_error = cudaFree(t1_dev);
  AssertCuda(cuda_error);
  cuda_error = cudaFree(t2_dev);
  AssertCuda(cuda_error);
}

int
main()
{
  initlog();

  init_cuda();

  test_gpu<0, 3, double>();
  test_gpu<1, 3, double>();
  test_gpu<2, 3, double>();
  test_gpu<0, 3, float>();
  test_gpu<1, 3, float>();
  test_gpu<2, 3, float>();
}
