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


// like _03 but for SymmetricTensor

#include "../tests.h"
#include <iomanip>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <typeinfo>
#include <complex>

#include <deal.II/base/template_constraints.h>
#include <deal.II/base/symmetric_tensor.h>


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
  check<SymmetricTensor<2,1,double>,double,SymmetricTensor<2,1,double> >();
  check<double,SymmetricTensor<2,1,double>,SymmetricTensor<2,1,double> >();

  check<SymmetricTensor<2,1,double>,float,SymmetricTensor<2,1,double> >();
  check<float,SymmetricTensor<2,1,double>,SymmetricTensor<2,1,double> >();

  check<SymmetricTensor<2,1,double>,std::complex<double>,SymmetricTensor<2,1,std::complex<double> > >();
  check<std::complex<double>,SymmetricTensor<2,1,double>,SymmetricTensor<2,1,std::complex<double> > >();

  check<SymmetricTensor<2,1,double>,std::complex<float>,SymmetricTensor<2,1,std::complex<double> > >();
  check<std::complex<float>,SymmetricTensor<2,1,double>,SymmetricTensor<2,1,std::complex<double> > >();

  check<SymmetricTensor<2,1,float>,std::complex<double>,SymmetricTensor<2,1,std::complex<double> > >();
  check<std::complex<double>,SymmetricTensor<2,1,float>,SymmetricTensor<2,1,std::complex<double> > >();

  check<SymmetricTensor<2,1,float>,std::complex<float>,SymmetricTensor<2,1,std::complex<float> > >();
  check<std::complex<float>,SymmetricTensor<2,1,float>,SymmetricTensor<2,1,std::complex<float> > >();
  
  deallog << "OK" << std::endl;
}
