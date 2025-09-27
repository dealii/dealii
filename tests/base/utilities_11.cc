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


// This tests the implementation of Utilities::to_string for different types.
// Note that the floating point number output might be depend on the system.

#include <deal.II/base/utilities.h>

#include "../tests.h"


void
test()
{
  const auto i = static_cast<unsigned long long int>(std::pow(2, 33));
  Assert(Utilities::to_string(i) == "8589934592", ExcInternalError());
  Assert(Utilities::to_string(i, 11) == "08589934592", ExcInternalError());

  const auto j = static_cast<unsigned long long int>(std::pow(2, 31));
  Assert(Utilities::to_string(j) == "2147483648", ExcInternalError());

  const int k = static_cast<int>(-std::pow(2, 30));
  Assert(Utilities::to_string(k) == "-1073741824", ExcInternalError());
  Assert(Utilities::to_string(k, 12) == "-01073741824", ExcInternalError());

  const auto l = static_cast<long long int>(-std::pow(2, 35));
  Assert(Utilities::to_string(l) == "-34359738368", ExcInternalError());
  Assert(Utilities::to_string(l, 13) == "-034359738368", ExcInternalError());

  float f(-3.14159265358979323846264338327950288419716939937510);
  Assert(Utilities::to_string(f) == "-3.14159274", ExcInternalError());
  Assert(Utilities::to_string(f, 13) == "-003.14159274", ExcInternalError());

  double d(-3.14159265358979323846264338327950288419716939937510);
  Assert(Utilities::to_string(d) == "-3.1415926535897931", ExcInternalError());
  Assert(Utilities::to_string(d, 20) == "-03.1415926535897931",
         ExcInternalError());

  long double ld(-3.14159265358979323846264338327950288419716939937510L);
  Assert(Utilities::to_string(ld) == "-3.14159265358979323851",
         ExcInternalError());
  Assert(Utilities::to_string(ld, 24) == "-03.14159265358979323851",
         ExcInternalError());

  double ed(-3.1415926535e-115);
  Assert(Utilities::to_string(ed) == "-3.1415926534999999e-115",
         ExcInternalError());
  Assert(Utilities::to_string(ed, 28) == "-00003.1415926534999999e-115",
         ExcInternalError());

  deallog << "Ok." << std::endl;
}



int
main()
{
  initlog();
  deallog.depth_console(0);

  test();
}
