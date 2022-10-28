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
miscellaneous_kernel(Number *check_1,
                     Number *check_2,
                     Number *check_3,
                     Number *check_4,
                     Number *check_5)
{
  // constructors
  typename Tensor<rank, dim, Number>::array_type array{};
  Tensor<rank, dim, Number>                      dummy_1(array);
  *check_1 = dummy_1.norm_square();
  Tensor<rank, dim, Number> dummy_2;
  *check_2                          = dummy_2.norm_square();
  Tensor<rank, dim, Number> dummy_3 = dummy_2;
  *check_3                          = dummy_3.norm_square();

  // access
  Tensor<rank + 1, dim, Number>   initializer_1;
  const Tensor<rank, dim, Number> dummy_5 = initializer_1[0];
  *check_4                                = dummy_5.norm_square();

  // assignment
  dummy_2  = dummy_3;
  *check_5 = dummy_2.norm_square();
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
__global__ void
init_kernel(Tensor<0, dim, Number> *t)
{
  if (threadIdx.x == 0)
    *t = 1.;
}

template <int dim, typename Number>
__global__ void
init_kernel(Tensor<1, dim, Number> *t)
{
  const unsigned int i = threadIdx.x;
  if (i < dim)
    (*t)[i] = i + 1.;
}

template <int dim, typename Number>
__global__ void
init_kernel(Tensor<2, dim, Number> *t)
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

  // Free memory
  cuda_error = cudaFree(t_dev);
  AssertCuda(cuda_error);
  cuda_error = cudaFree(t1_dev);
  AssertCuda(cuda_error);
  cuda_error = cudaFree(t2_dev);
  AssertCuda(cuda_error);

  // Miscellaneous
  {
    Number *check_1;
    Number *check_2;
    Number *check_3;
    Number *check_4;
    Number *check_5;

    cuda_error = cudaMalloc(&check_1, sizeof(Number));
    AssertCuda(cuda_error);
    cuda_error = cudaMalloc(&check_2, sizeof(Number));
    AssertCuda(cuda_error);
    cuda_error = cudaMalloc(&check_3, sizeof(Number));
    AssertCuda(cuda_error);
    cuda_error = cudaMalloc(&check_4, sizeof(Number));
    AssertCuda(cuda_error);
    cuda_error = cudaMalloc(&check_5, sizeof(Number));
    AssertCuda(cuda_error);

    miscellaneous_kernel<rank, dim, Number>
      <<<1, 1>>>(check_1, check_2, check_3, check_4, check_5);

    Number check_1_host, check_2_host, check_3_host, check_4_host, check_5_host;

    cuda_error = cudaMemcpy(&check_1_host,
                            check_1,
                            sizeof(Number),
                            cudaMemcpyDeviceToHost);
    AssertCuda(cuda_error);
    cuda_error = cudaMemcpy(&check_2_host,
                            check_2,
                            sizeof(Number),
                            cudaMemcpyDeviceToHost);
    AssertCuda(cuda_error);
    cuda_error = cudaMemcpy(&check_3_host,
                            check_3,
                            sizeof(Number),
                            cudaMemcpyDeviceToHost);
    AssertCuda(cuda_error);
    cuda_error = cudaMemcpy(&check_4_host,
                            check_4,
                            sizeof(Number),
                            cudaMemcpyDeviceToHost);
    AssertCuda(cuda_error);
    cuda_error = cudaMemcpy(&check_5_host,
                            check_5,
                            sizeof(Number),
                            cudaMemcpyDeviceToHost);
    AssertCuda(cuda_error);

    AssertThrow(std::abs(check_1_host) < tolerance, ExcInternalError());
    AssertThrow(std::abs(check_2_host) < tolerance, ExcInternalError());
    AssertThrow(std::abs(check_3_host) < tolerance, ExcInternalError());
    AssertThrow(std::abs(check_4_host) < tolerance, ExcInternalError());
    AssertThrow(std::abs(check_5_host) < tolerance, ExcInternalError());

    cuda_error = cudaFree(check_1);
    AssertCuda(cuda_error);
    cuda_error = cudaFree(check_2);
    AssertCuda(cuda_error);
    cuda_error = cudaFree(check_3);
    AssertCuda(cuda_error);
    cuda_error = cudaFree(check_4);
    AssertCuda(cuda_error);
    cuda_error = cudaFree(check_5);
    AssertCuda(cuda_error);
  }
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
