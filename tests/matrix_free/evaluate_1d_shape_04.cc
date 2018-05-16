// ---------------------------------------------------------------------
//
// Copyright (C) 2018 by the deal.II authors
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



// check the correctness of the 1d evaluation functions used in FEEvaluation,
// path evaluate_general, when using same array for in and out

#include "../tests.h"
#include <iostream>

#include <deal.II/matrix_free/tensor_product_kernels.h>


template <int M, int N, int type>
void
test()
{
  deallog << "Test " << M << " x " << N << std::endl;
  AlignedVector<double> shape(M*N);
  for (unsigned int i=0; i<M; ++i)
    for (unsigned int j=0; j<N; ++j)
      shape[i*N+j] = -1. + 2. * random_value<double>();

  double x[N+M], x_ref[N], y_ref[M];
  for (unsigned int i=0; i<N; ++i)
    x[i] = random_value<double>();

  // compute reference
  for (unsigned int i=0; i<M; ++i)
    {
      y_ref[i] = 0.;
      for (unsigned int j=0; j<N; ++j)
        y_ref[i] += shape[i*N+j] * x[j];
    }

  // apply function for tensor product
  internal::EvaluatorTensorProduct<internal::evaluate_general,1,M,N,double>
  evaluator(shape, shape, shape);
  if (type == 0)
    evaluator.template values<0,false,false> (x,x);
  if (type == 1)
    evaluator.template gradients<0,false,false> (x,x);
  if (type == 2)
    evaluator.template hessians<0,false,false> (x,x);

  deallog << "Errors no transpose: ";
  for (unsigned int i=0; i<M; ++i)
    deallog << x[i] - y_ref[i] << " ";
  deallog << std::endl;


  for (unsigned int i=0; i<M; ++i)
    x[i] = random_value<double>();

  // compute reference
  for (unsigned int i=0; i<N; ++i)
    {
      x_ref[i] = 0.;
      for (unsigned int j=0; j<M; ++j)
        x_ref[i] += shape[j*N+i] * x[j];
    }

  // apply function for tensor product
  if (type == 0)
    evaluator.template values<0,true,false> (x,x);
  if (type == 1)
    evaluator.template gradients<0,true,false> (x,x);
  if (type == 2)
    evaluator.template hessians<0,true,false> (x,x);

  deallog << "Errors transpose:    ";
  for (unsigned int i=0; i<N; ++i)
    deallog << x[i] - x_ref[i] << " ";
  deallog << std::endl;
}

int
main ()
{
  initlog();

  deallog.push("values");
  test<4,4,0>();
  test<3,3,0>();
  test<4,3,0>();
  test<3,4,0>();
  test<3,5,0>();
  deallog.pop();

  deallog.push("gradients");
  test<4,4,1>();
  test<3,3,1>();
  test<4,3,1>();
  test<3,4,1>();
  test<3,5,1>();
  deallog.pop();

  deallog.push("hessians");
  test<4,4,2>();
  test<3,3,2>();
  test<4,3,2>();
  test<3,4,2>();
  test<3,5,2>();
  deallog.pop();

  return 0;
}
