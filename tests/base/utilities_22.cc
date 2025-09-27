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


// test Utilities::pow where we should get an integer overflow and thus a
// compile error

#include <deal.II/base/utilities.h>

#include "../tests.h"



int
main()
{
  initlog();

  constexpr int a = Utilities::pow(216, 4);
  deallog << a << std::endl;
}
