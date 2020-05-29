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
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------



// this function tests the correctness of the 2d evaluation functions used in
// CUDAWrappers::FEEvaluation. These functions are marked 'internal' but it is
// much easier to check their correctness directly rather than from the results
// in dependent functions

#include <deal.II/lac/cuda_vector.h>
#include <deal.II/lac/read_write_vector.h>

#include <deal.II/matrix_free/cuda_fe_evaluation.h>

#include <fstream>
#include <iostream>

#include "../tests.h"

namespace CUDA = LinearAlgebra::CUDAWrappers;

template <int M, int N, int type, bool add, bool dof_to_quad>
__global__ void
evaluate_tensor_product(double *dst, double *src)
{
  CUDAWrappers::internal::EvaluatorTensorProduct<
    CUDAWrappers::internal::evaluate_general,
    2,
    M - 1,
    N,
    double>
    evaluator(0);

  if (type == 0)
    {
      evaluator.template values<0, dof_to_quad, false, false>(
        CUDAWrappers::internal::get_global_shape_values<double>(0), src, src);
      __syncthreads();
      evaluator.template values<1, dof_to_quad, add, false>(
        CUDAWrappers::internal::get_global_shape_values<double>(0), src, dst);
    }
  if (type == 1)
    {
      evaluator.template gradients<0, dof_to_quad, false, false>(
        CUDAWrappers::internal::get_global_shape_values<double>(0), src, src);
      __syncthreads();
      evaluator.template gradients<1, dof_to_quad, add, false>(
        CUDAWrappers::internal::get_global_shape_values<double>(0), src, dst);
    }
}

template <int M, int N, int type, bool add>
void
test()
{
  deallog << "Test " << M << " x " << N << std::endl;
  LinearAlgebra::ReadWriteVector<double> shape_host(M * N);
  for (unsigned int i = 0; i < (M + 1) / 2; ++i)
    for (unsigned int j = 0; j < N; ++j)
      {
        shape_host[i * N + j] =
          -1. + 2. * static_cast<double>(Testing::rand()) / RAND_MAX;
        if (type == 1)
          shape_host[(M - 1 - i) * N + N - 1 - j] = -shape_host[i * N + j];
        else
          shape_host[(M - 1 - i) * N + N - 1 - j] = shape_host[i * N + j];
      }
  if (type == 0 && M % 2 == 1 && N % 2 == 1)
    {
      for (unsigned int i = 0; i < M; ++i)
        shape_host[i * N + N / 2] = 0.;
      shape_host[M / 2 * N + N / 2] = 1.;
    }
  if (type == 1 && M % 2 == 1 && N % 2 == 1)
    shape_host[M / 2 * N + N / 2] = 0.;

  constexpr int                          M_2d = M * M;
  constexpr int                          N_2d = N * N;
  LinearAlgebra::ReadWriteVector<double> x_host(N_2d), x_ref(N_2d),
    y_host(M_2d), y_ref(M_2d);
  for (unsigned int i = 0; i < N_2d; ++i)
    x_host[i] = static_cast<double>(Testing::rand()) / RAND_MAX;

  FullMatrix<double> shape_2d(M_2d, N_2d);
  for (unsigned int i = 0; i < M; ++i)
    {
      for (unsigned int j = 0; j < N; ++j)
        {
          const double shape_val = shape_host[i * N + j];
          for (unsigned int m = 0; m < M; ++m)
            for (unsigned int n = 0; n < N; ++n)
              shape_2d(i * M + m, j * N + n) =
                shape_val * shape_host[m * N + n];
        }
    }

  // Compute reference
  for (unsigned int i = 0; i < M_2d; ++i)
    {
      y_host[i] = 1.;
      y_ref[i]  = add ? y_host[i] : 0.;
      for (unsigned int j = 0; j < N_2d; ++j)
        y_ref[i] += shape_2d(i, j) * x_host[j];
    }

  // Copy data to the GPU.
  CUDA::Vector<double> x_dev(N_2d), y_dev(M_2d);
  x_dev.import(x_host, VectorOperation::insert);
  y_dev.import(y_host, VectorOperation::insert);

  unsigned int size_shape_values = M * N * sizeof(double);

  cudaError_t cuda_error =
    cudaMemcpyToSymbol(CUDAWrappers::internal::get_global_shape_values<double>(
                         0),
                       shape_host.begin(),
                       size_shape_values,
                       0,
                       cudaMemcpyHostToDevice);
  AssertCuda(cuda_error);

  cuda_error = cudaMemcpyToSymbol(
    CUDAWrappers::internal::get_global_shape_gradients<double>(0),
    shape_host.begin(),
    size_shape_values,
    0,
    cudaMemcpyHostToDevice);
  AssertCuda(cuda_error);

  // Launch the kernel
  dim3 block_dim(M, N);
  evaluate_tensor_product<M, N, type, add, false>
    <<<1, block_dim>>>(y_dev.get_values(), x_dev.get_values());

  // Check the results on the host
  y_host.import(y_dev, VectorOperation::insert);
  deallog << "Errors no transpose: ";

  for (unsigned int i = 0; i < M_2d; ++i)
    deallog << y_host[i] - y_ref[i] << " ";
  deallog << std::endl;

  for (unsigned int i = 0; i < M_2d; ++i)
    y_host[i] = static_cast<double>(Testing::rand()) / RAND_MAX;

  // Copy y_host to the device
  y_dev.import(y_host, VectorOperation::insert);

  // Compute reference
  for (unsigned int i = 0; i < N_2d; ++i)
    {
      x_host[i] = 2.;
      x_ref[i]  = add ? x_host[i] : 0.;
      for (unsigned int j = 0; j < M_2d; ++j)
        x_ref[i] += shape_2d(j, i) * y_host[j];
    }

  // Copy x_host to the device
  x_dev.import(x_host, VectorOperation::insert);

  // Launch the kernel
  evaluate_tensor_product<M, N, type, add, true>
    <<<1, block_dim>>>(x_dev.get_values(), y_dev.get_values());

  // Check the results on the host
  x_host.import(x_dev, VectorOperation::insert);
  deallog << "Errors transpose:    ";
  for (unsigned int i = 0; i < N_2d; ++i)
    deallog << x_host[i] - x_ref[i] << " ";
  deallog << std::endl;
}

int
main()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);

  init_cuda();

  deallog.push("values");
  test<4, 4, 0, false>();
  test<3, 3, 0, false>();
  deallog.pop();

  deallog.push("gradients");
  test<4, 4, 1, false>();
  test<3, 3, 1, false>();
  deallog.pop();

  deallog.push("add");

  deallog.push("values");
  test<4, 4, 0, true>();
  test<3, 3, 0, true>();
  deallog.pop();

  deallog.push("gradients");
  test<4, 4, 1, true>();
  test<3, 3, 1, true>();
  deallog.pop();

  deallog.pop();

  return 0;
}
