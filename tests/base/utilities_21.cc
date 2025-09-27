// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test Utilities::pow

#include <deal.II/base/utilities.h>

#include "../tests.h"



int
main()
{
  initlog();

  constexpr int a = Utilities::pow(100, 2);
  constexpr int b = Utilities::pow(200, 2);
  constexpr int c = Utilities::pow(215, 2);
  constexpr int d = Utilities::pow(216, 2);
  constexpr int e = Utilities::pow(500, 2);
  constexpr int f = Utilities::pow(600, 2);
  constexpr int g = Utilities::pow(700, 2);
  constexpr int h = Utilities::pow(800, 2);
  constexpr int i = Utilities::pow(900, 2);
  constexpr int j = Utilities::pow(215, 3);
  constexpr int k = Utilities::pow(216, 3);
  constexpr int l = Utilities::pow(1285, 3);
  deallog << a + b + c + d + e + f + g + h + i + j + k + l << std::endl;

  return 0;
}
