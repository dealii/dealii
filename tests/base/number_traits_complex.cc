// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2017 by the deal.II authors
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
  std::ofstream logfile("output");
  deallog << std::setprecision(3);
  deallog.attach(logfile);

  check(std::complex<float>(1.5, 2.5));
  check(std::complex<float>(-1.5, -2.5));

  check(std::complex<double>(1.5, 2.5));
  check(std::complex<double>(-1.5, -2.5));

  check(std::complex<long double>(1.5, 2.5));
  check(std::complex<long double>(-1.5, -2.5));

  return 0;
}
