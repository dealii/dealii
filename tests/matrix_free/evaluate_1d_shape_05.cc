// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check the correctness of the 1d evaluation functions used in FEEvaluation,
// path evaluate_evenodd, when using same array for in and out

#include <deal.II/matrix_free/tensor_product_kernels.h>

#include <iostream>

#include "../tests.h"


template <int M, int N, int type>
void
test()
{
  deallog << "Test " << M << " x " << N << std::endl;
  double shape[M][N];
  for (unsigned int i = 0; i < (M + 1) / 2; ++i)
    for (unsigned int j = 0; j < N; ++j)
      {
        shape[i][j] = -1. + 2. * random_value<double>();
        if (type == 1)
          shape[M - 1 - i][N - 1 - j] = -shape[i][j];
        else
          shape[M - 1 - i][N - 1 - j] = shape[i][j];
      }
  if (type == 0 && M % 2 == 1 && N % 2 == 1)
    {
      for (unsigned int i = 0; i < M; ++i)
        shape[i][N / 2] = 0.;
      shape[M / 2][N / 2] = 0.9;
    }
  if (type == 1 && M % 2 == 1 && N % 2 == 1)
    shape[M / 2][N / 2] = 0.;


  // create symmetrized shape array exactly as expected by the evenodd
  // function
  AlignedVector<double> shape_sym(M * ((N + 1) / 2));
  for (unsigned int i = 0; i < M / 2; ++i)
    for (unsigned int q = 0; q < (N + 1) / 2; ++q)
      {
        shape_sym[i * ((N + 1) / 2) + q] =
          0.5 * (shape[i][q] + shape[i][N - 1 - q]);
        shape_sym[(M - 1 - i) * ((N + 1) / 2) + q] =
          0.5 * (shape[i][q] - shape[i][N - 1 - q]);
      }
  if (M % 2 == 1)
    for (unsigned int q = 0; q < (N + 1) / 2; ++q)
      shape_sym[(M - 1) / 2 * ((N + 1) / 2) + q] = shape[(M - 1) / 2][q];

  double x[N + M], x_ref[N], y_ref[M];
  for (unsigned int i = 0; i < N; ++i)
    x[i] = random_value<double>();

  // compute reference
  for (unsigned int i = 0; i < M; ++i)
    {
      y_ref[i] = 0.;
      for (unsigned int j = 0; j < N; ++j)
        y_ref[i] += shape[i][j] * x[j];
    }

  // apply function for tensor product
  internal::EvaluatorTensorProduct<internal::evaluate_evenodd, 1, M, N, double>
    evaluator(shape_sym, shape_sym, shape_sym);
  if (type == 0)
    evaluator.template values<0, false, false>(x, x);
  if (type == 1)
    evaluator.template gradients<0, false, false>(x, x);
  if (type == 2)
    evaluator.template hessians<0, false, false>(x, x);

  deallog << "Errors no transpose: ";
  for (unsigned int i = 0; i < M; ++i)
    deallog << x[i] - y_ref[i] << ' ';
  deallog << std::endl;


  for (unsigned int i = 0; i < M; ++i)
    x[i] = random_value<double>();

  // compute reference
  for (unsigned int i = 0; i < N; ++i)
    {
      x_ref[i] = 0.;
      for (unsigned int j = 0; j < M; ++j)
        x_ref[i] += shape[j][i] * x[j];
    }

  // apply function for tensor product
  if (type == 0)
    evaluator.template values<0, true, false>(x, x);
  if (type == 1)
    evaluator.template gradients<0, true, false>(x, x);
  if (type == 2)
    evaluator.template hessians<0, true, false>(x, x);

  deallog << "Errors transpose:    ";
  for (unsigned int i = 0; i < N; ++i)
    deallog << x[i] - x_ref[i] << ' ';
  deallog << std::endl;
}

int
main()
{
  initlog();

  deallog.push("values");
  test<4, 4, 0>();
  test<3, 3, 0>();
  test<4, 3, 0>();
  test<3, 4, 0>();
  test<3, 5, 0>();
  deallog.pop();

  deallog.push("gradients");
  test<4, 4, 1>();
  test<3, 3, 1>();
  test<4, 3, 1>();
  test<3, 4, 1>();
  test<3, 5, 1>();
  deallog.pop();

  deallog.push("hessians");
  test<4, 4, 2>();
  test<3, 3, 2>();
  test<4, 3, 2>();
  test<3, 4, 2>();
  test<3, 5, 2>();
  deallog.pop();

  return 0;
}
