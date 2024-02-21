// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test the implementation of theh constexpr min/max functions
// in case deal.II is compiled with C++14 support, we only test
// the STL.


#include <algorithm>

#include "../tests.h"

constexpr bool
comp(const int &a, const int &b)
{
  return b < a;
}


int
main()
{
  initlog();

  constexpr int max_1 = std::max(0, 1);
  deallog << max_1 << std::endl;
  constexpr int max_2 = std::max(3, 2, comp);
  deallog << max_2 << std::endl;

  constexpr int min_1 = std::min(1, 2);
  deallog << min_1 << std::endl;
  constexpr int min_2 = std::min(1, 2, comp);
  deallog << min_2 << std::endl;

  deallog << "OK" << std::endl;
}
