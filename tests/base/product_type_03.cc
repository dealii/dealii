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


// test that the product between a tensor of rank 2 and a scalar
// results a result type as expected

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
  Assert (typeid(T() * U()) == typeid(CompareType),
	  ExcInternalError());
  Assert (typeid(T() * U()) == typeid(CompareType),
	  ExcInternalError());
}


int main()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  // check product of scalars
  check<Tensor<2,1,double>,double,Tensor<2,1,double> >();
  check<double,Tensor<2,1,double>,Tensor<2,1,double> >();

  check<Tensor<2,1,double>,float,Tensor<2,1,double> >();
  check<float,Tensor<2,1,double>,Tensor<2,1,double> >();

  check<Tensor<2,1,double>,std::complex<double>,Tensor<2,1,std::complex<double> > >();
  check<std::complex<double>,Tensor<2,1,double>,Tensor<2,1,std::complex<double> > >();

  check<Tensor<2,1,double>,std::complex<float>,Tensor<2,1,std::complex<double> > >();
  check<std::complex<float>,Tensor<2,1,double>,Tensor<2,1,std::complex<double> > >();

  check<Tensor<2,1,float>,std::complex<double>,Tensor<2,1,std::complex<double> > >();
  check<std::complex<double>,Tensor<2,1,float>,Tensor<2,1,std::complex<double> > >();

  check<Tensor<2,1,float>,std::complex<float>,Tensor<2,1,std::complex<float> > >();
  check<std::complex<float>,Tensor<2,1,float>,Tensor<2,1,std::complex<float> > >();
  
  deallog << "OK" << std::endl;
}
