// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2017 by the deal.II authors
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

// Test for mixed Number type operations of Tensor<2, dim, Number> and
// Tensor<3, dim. Number> (addition, subtraction, contraction).

#include "../tests.h"
#include <deal.II/base/tensor.h>

#include<complex>


int
main ()
{
  initlog();

  float f_scalar = 10.;
  double d_scalar = 10.;
  std::complex<double> c_scalar = 10.;

  Tensor<1,2,float>                 f1;
  f1[0] = 1.;
  f1[1] = 2.;
  Tensor<1,2,double>                d1;
  d1[0] = 4.;
  d1[1] = 8.;
  Tensor<1,2,std::complex<double> > c1;
  c1[0] = 16.;
  c1[1] = 32.;

  deallog << f1 + d1 << std::endl;
  deallog << f1 - d1 << std::endl;
  deallog << f1 *d1 << std::endl;

  deallog << d1 + c1 << std::endl;
  deallog << d1 - c1 << std::endl;
  deallog << d1 *c1 << std::endl;

  deallog << (d1 * f_scalar) << std::endl;
  deallog << (d1 / f_scalar) << std::endl;
  deallog << (f_scalar * d1) << std::endl;
  deallog << (d1 *= f_scalar) << std::endl;
  deallog << (d1 /= f_scalar) << std::endl;

  deallog << (c1 * d_scalar) << std::endl;
  deallog << (c1 / d_scalar) << std::endl;
  deallog << (d_scalar * c1) << std::endl;
  deallog << (c1 *= d_scalar) << std::endl;
  deallog << (c1 /= d_scalar) << std::endl;

  Tensor<2,2,float> f2;
  f2[0] = f1;
  f2[1] = 2. * f1;

  Tensor<2,2,double> d2;
  d2[0] = d1;
  d2[1] = 2. * d1;

  Tensor<2,2,std::complex<double> > c2;
  c2[0] = c1;
  c2[1] = 2. * c1;

  deallog << f2 + d2 << std::endl;
  deallog << f2 - d2 << std::endl;
  deallog << f2 *d2 << std::endl;

  deallog << d2 + c2 << std::endl;
  deallog << d2 - c2 << std::endl;
  deallog << d2 *c2 << std::endl;

  deallog << (d2 * f_scalar) << std::endl;
  deallog << (d2 / f_scalar) << std::endl;
  deallog << (f_scalar * d2) << std::endl;
  deallog << (d2 *= f_scalar) << std::endl;
  deallog << (d2 /= f_scalar) << std::endl;

  deallog << (c2 * d_scalar) << std::endl;
  deallog << (c2 / d_scalar) << std::endl;
  deallog << (d_scalar * c2) << std::endl;
  deallog << (c2 *= d_scalar) << std::endl;
  deallog << (c2 /= d_scalar) << std::endl;

  Tensor<3,2,float> f3;
  f3[0] = f2;
  f3[1] = 2. * f2;

  Tensor<3,2,double> d3;
  d3[0] = d2;
  d3[1] = 2. * d2;

  Tensor<3,2,std::complex<double> > c3;
  c3[0] = c2;
  c3[1] = 2. * c2;

  deallog << f3 + d3 << std::endl;
  deallog << f3 - d3 << std::endl;
  deallog << f3 *d3 << std::endl;

  deallog << d3 + c3 << std::endl;
  deallog << d3 - c3 << std::endl;
  deallog << d3 *c3 << std::endl;

  deallog << (d3 * f_scalar) << std::endl;
  deallog << (d3 / f_scalar) << std::endl;
  deallog << (f_scalar * d3) << std::endl;
  deallog << (d3 *= f_scalar) << std::endl;
  deallog << (d3 /= f_scalar) << std::endl;

  deallog << (c3 * d_scalar) << std::endl;
  deallog << (c3 / d_scalar) << std::endl;
  deallog << (d_scalar * c3) << std::endl;
  deallog << (c3 *= d_scalar) << std::endl;
  deallog << (c3 /= d_scalar) << std::endl;

  deallog << "OK!" << std::endl;
}
