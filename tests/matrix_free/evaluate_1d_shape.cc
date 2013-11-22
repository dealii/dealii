// ---------------------------------------------------------------------
// $Id$
//
// Copyright (C) 2013 by the deal.II authors
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



// this function tests the correctness of the 1d evaluation functions used in
// FEEvaluation. These functions are marked 'internal' but it is much easier
// to check their correctness directly rather than from the results in
// dependent functions

#include "../tests.h"
#include <iostream>
#include <fstream>

#include <deal.II/matrix_free/fe_evaluation.h>


template <int M, int N, int type, bool add>
void test()
{
  deallog << "Test " << M << " x " << N << std::endl;
  double shape[M][N];
  for (unsigned int i=0; i<(M+1)/2; ++i)
    for (unsigned int j=0; j<N; ++j)
      {
        shape[i][j] = -1. + 2. * (double)Testing::rand()/RAND_MAX;
        if (type == 1)
          shape[M-1-i][N-1-j] = -shape[i][j];
        else
          shape[M-1-i][N-1-j] = shape[i][j];
      }
  if (type == 0 && M%2 == 1 && N%2 == 1)
    {
      for (unsigned int i=0; i<M; ++i)
        shape[i][N/2] = 0.;
      shape[M/2][N/2] = 1;
    }
  if (type == 1 && M%2 == 1 && N%2 == 1)
    shape[M/2][N/2] = 0.;

  double x[N], x_ref[N], y[M], y_ref[M];
  for (unsigned int i=0; i<N; ++i)
    x[i] = (double)Testing::rand()/RAND_MAX;

  // compute reference
  for (unsigned int i=0; i<M; ++i)
    {
      y[i] = 1.;
      y_ref[i] = add ? y[i] : 0.;
      for (unsigned int j=0; j<N; ++j)
        y_ref[i] += shape[i][j] * x[j];
    }

  // apply function for tensor product
  if (type == 0)
    internal::apply_tensor_product_values<1,M-1,N,double,0,false,add>
      (&shape[0][0],x,y);
  if (type == 1)
    internal::apply_tensor_product_gradients<1,M-1,N,double,0,false,add>
      (&shape[0][0],x,y);
  if (type == 2)
    internal::apply_tensor_product_hessians<1,M-1,N,double,0,false,add>
      (&shape[0][0],x,y);

  deallog << "Errors no transpose: ";
  for (unsigned int i=0; i<M; ++i)
    deallog << y[i] - y_ref[i] << " ";
  deallog << std::endl;


  for (unsigned int i=0; i<M; ++i)
    y[i] = (double)Testing::rand()/RAND_MAX;

  // compute reference
  for (unsigned int i=0; i<N; ++i)
    {
      x[i] = 2.;
      x_ref[i] = add ? x[i] : 0.;
      for (unsigned int j=0; j<M; ++j)
        x_ref[i] += shape[j][i] * y[j];
    }

  // apply function for tensor product
  if (type == 0)
    internal::apply_tensor_product_values<1,M-1,N,double,0,true,add>
      (&shape[0][0],y,x);
  if (type == 1)
    internal::apply_tensor_product_gradients<1,M-1,N,double,0,true,add>
      (&shape[0][0],y,x);
  if (type == 2)
    internal::apply_tensor_product_hessians<1,M-1,N,double,0,true,add>
      (&shape[0][0],y,x);

  deallog << "Errors transpose:    ";
  for (unsigned int i=0; i<N; ++i)
    deallog << x[i] - x_ref[i] << " ";
  deallog << std::endl;
}

int main ()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1e-14);

  deallog.push("values");
  test<4,4,0,false>();
  test<3,3,0,false>();
  test<4,3,0,false>();
  test<3,4,0,false>();
  test<3,5,0,false>();
  deallog.pop();

  deallog.push("gradients");
  test<4,4,1,false>();
  test<3,3,1,false>();
  test<4,3,1,false>();
  test<3,4,1,false>();
  test<3,5,1,false>();
  deallog.pop();

  deallog.push("hessians");
  test<4,4,2,false>();
  test<3,3,2,false>();
  test<4,3,2,false>();
  test<3,4,2,false>();
  test<3,5,2,false>();
  deallog.pop();

  deallog.push("add");

  deallog.push("values");
  test<4,4,0,true>();
  test<3,3,0,true>();
  test<4,3,0,true>();
  test<3,4,0,true>();
  test<3,5,0,true>();
  deallog.pop();

  deallog.push("gradients");
  test<4,4,1,true>();
  test<3,3,1,true>();
  test<4,3,1,true>();
  test<3,4,1,true>();
  test<3,5,1,true>();
  deallog.pop();

  deallog.push("hessians");
  test<4,4,2,true>();
  test<3,3,2,true>();
  test<4,3,2,true>();
  test<3,4,2,true>();
  test<3,5,2,true>();
  deallog.pop();

  deallog.pop();

  return 0;
}
    
