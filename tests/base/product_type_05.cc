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


// test that ProductType<double,double> can be resolved. same for
// ProductType<std::complex<double>,std::complex<double> >

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
  Assert (typeid(typename ProductType<T,U>::type) == typeid(CompareType),
	  ExcInternalError());
  Assert (typeid(typename ProductType<T,U>::type) == typeid(T() * U()),
	  ExcInternalError());
}


int main()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  check<double,double,double>();
  deallog << typename ProductType<double,double>::type(2.345)*typename ProductType<double,double>::type(3.456)
	  << ' '
	  << typename ProductType<double,double>::type(2.345*3.456)
	  << std::endl;
  
  check<std::complex<double>,std::complex<double>,std::complex<double> >();
  deallog << (typename ProductType<std::complex<double>,std::complex<double> >::type(2.345, 1.23) *
 	      typename ProductType<std::complex<double>,std::complex<double> >::type(3.456, 2.45))
	  << ' '
  	  << (typename ProductType<std::complex<double>,std::complex<double> >::type
  	      (std::complex<double>(2.345, 1.23) *
  	       std::complex<double>(3.456, 2.45)))
  	  << std::endl;
  
  deallog << "OK" << std::endl;
}
