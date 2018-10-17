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
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------


// This tests the implementation of Utilities::to_string for different types.
// Note that the floating point number output might be depend on the system.

#include <deal.II/base/utilities.h>

#include "../tests.h"


void
test()
{
  unsigned long long int i = std::pow(2, 33);
  DEAL_II_Assert(Utilities::to_string(i) == "8589934592", ExcInternalError());
  DEAL_II_Assert(Utilities::to_string(i, 11) == "08589934592",
                 ExcInternalError());

  unsigned long long j = std::pow(2, 31);
  DEAL_II_Assert(Utilities::to_string(j) == "2147483648", ExcInternalError());

  int k = -std::pow(2, 30);
  DEAL_II_Assert(Utilities::to_string(k) == "-1073741824", ExcInternalError());
  DEAL_II_Assert(Utilities::to_string(k, 12) == "-01073741824",
                 ExcInternalError());

  long long int l = -std::pow(2, 35);
  DEAL_II_Assert(Utilities::to_string(l) == "-34359738368", ExcInternalError());
  DEAL_II_Assert(Utilities::to_string(l, 13) == "-034359738368",
                 ExcInternalError());

  float f(-3.14159265358979323846264338327950288419716939937510);
  DEAL_II_Assert(Utilities::to_string(f) == "-3.14159274", ExcInternalError());
  DEAL_II_Assert(Utilities::to_string(f, 13) == "-003.14159274",
                 ExcInternalError());

  double d(-3.14159265358979323846264338327950288419716939937510);
  DEAL_II_Assert(Utilities::to_string(d) == "-3.1415926535897931",
                 ExcInternalError());
  DEAL_II_Assert(Utilities::to_string(d, 20) == "-03.1415926535897931",
                 ExcInternalError());

  long double ld(-3.14159265358979323846264338327950288419716939937510L);
  DEAL_II_Assert(Utilities::to_string(ld) == "-3.14159265358979323851",
                 ExcInternalError());
  DEAL_II_Assert(Utilities::to_string(ld, 24) == "-03.14159265358979323851",
                 ExcInternalError());

  double ed(-3.1415926535e-115);
  DEAL_II_Assert(Utilities::to_string(ed) == "-3.1415926534999999e-115",
                 ExcInternalError());
  DEAL_II_Assert(Utilities::to_string(ed, 28) == "-00003.1415926534999999e-115",
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
