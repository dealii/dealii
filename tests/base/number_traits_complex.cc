// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2007 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// check numbers::NumberTraits for real data types

#include <limits>
#include <typeinfo>

#include "../tests.h"

// replace type names found on MAC OS
std::string
cleanup_type(const std::string &in)
{
  std::string ret = in;
  ret =
    Utilities::replace_in_string(ret, "NSt3__17complexIfEE", "St7complexIfE");
  ret =
    Utilities::replace_in_string(ret, "NSt3__17complexIdEE", "St7complexIdE");
  ret =
    Utilities::replace_in_string(ret, "NSt3__17complexIeEE", "St7complexIeE");
  return ret;
}


template <typename number>
void
check(const number &x)
{
  deallog << "typeid(x).name() = " << cleanup_type(typeid(x).name())
          << std::endl;

  deallog << "typeid(NumberTraits<number>::real_type).name() = "
          << typeid(typename numbers::NumberTraits<number>::real_type).name()
          << std::endl;

  deallog << numbers::NumberTraits<number>::conjugate(x) << std::endl;

  deallog << numbers::NumberTraits<number>::abs_square(x) << std::endl;

  deallog << numbers::NumberTraits<number>::abs(x) << std::endl;
}



int
main()
{
  initlog();
  deallog << std::setprecision(3);

  check(std::complex<float>(1.5, 2.5));
  check(std::complex<float>(-1.5, -2.5));

  check(std::complex<double>(1.5, 2.5));
  check(std::complex<double>(-1.5, -2.5));

  check(std::complex<long double>(1.5, 2.5));
  check(std::complex<long double>(-1.5, -2.5));

  return 0;
}
