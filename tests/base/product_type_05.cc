// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


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
  deallog
    << (ProductType<std::complex<double>, std::complex<double>>::type(2.345,
                                                                      1.23) *
        ProductType<std::complex<double>, std::complex<double>>::type(3.456,
                                                                      2.45))
    << ' '
    << (ProductType<std::complex<double>, std::complex<double>>::type(
         std::complex<double>(2.345, 1.23) * std::complex<double>(3.456, 2.45)))
    << std::endl;

  deallog << "OK" << std::endl;
}
