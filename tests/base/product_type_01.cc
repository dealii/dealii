// ---------------------------------------------------------------------
//
// Copyright (C) 2015 by the deal.II authors
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


// test ProductType

#include "../tests.h"
#include <iomanip>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <typeinfo>
#include <complex>

#include <deal.II/base/template_constraints.h>
#include <deal.II/base/tensor.h>


template <typename T, typename U, typename CompareType>
void check()
{
  AssertThrow (typeid(typename ProductType<T,U>::type) == typeid(CompareType),
               ExcInternalError());
  AssertThrow (typeid(typename ProductType<T,U>::type) == typeid(T() * U()),
               ExcInternalError());
}


int main()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  // check product of scalars
  check<double,double,double>();
  check<float,double,double>();
  check<double,float,double>();

  check<int,int,int>();
  check<int,double,double>();
  check<double,int,double>();

  // check product with Tensor<1,dim>
  check<Tensor<1,2,double>,double,Tensor<1,2,double> >();
  check<Tensor<1,2,float>,double,Tensor<1,2,double> >();
  check<double,Tensor<1,2,float>,Tensor<1,2,double> >();

  // check product with Tensor<2,dim> (which is a different class)
  check<Tensor<2,2,double>,double,Tensor<2,2,double> >();
  check<Tensor<2,2,float>,double,Tensor<2,2,double> >();
  check<double,Tensor<2,2,float>,Tensor<2,2,double> >();

  // check product with std::complex. rather annoyingly, there is no
  // product between std::complex<double> and float, or the other way
  // around, so stay within the same type system
  check<std::complex<double>,double,std::complex<double> >();
  check<std::complex<float>,float,std::complex<float> >();
  check<Tensor<1,2>,std::complex<double>,Tensor<1,2,std::complex<double> > >();
  check<std::complex<double>,Tensor<1,2>,Tensor<1,2,std::complex<double> > >();

  deallog << "OK" << std::endl;
}
