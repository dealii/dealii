// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2017 by the deal.II authors
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

#include <deal.II/base/template_constraints.h>
#include <deal.II/base/tensor.h>

#include <complex>
#include <typeinfo>

#include "../tests.h"



template <typename T, typename U, typename CompareType>
void
check()
{
  AssertThrow(typeid(typename ProductType<T, U>::type) == typeid(CompareType),
              ExcInternalError());
  AssertThrow(typeid(typename ProductType<T, U>::type) == typeid(T() * U()),
              ExcInternalError());
}


int
main()
{
  initlog();

  check<double, double, double>();
  deallog << ProductType<double, double>::type(2.345) *
               ProductType<double, double>::type(3.456)
          << ' ' << ProductType<double, double>::type(2.345 * 3.456)
          << std::endl;

  check<std::complex<double>, std::complex<double>, std::complex<double>>();
  deallog << (ProductType<std::complex<double>, std::complex<double>>::type(
                2.345, 1.23) *
              ProductType<std::complex<double>, std::complex<double>>::type(
                3.456, 2.45))
          << ' '
          << (ProductType<std::complex<double>, std::complex<double>>::type(
               std::complex<double>(2.345, 1.23) *
               std::complex<double>(3.456, 2.45)))
          << std::endl;

  deallog << "OK" << std::endl;
}
