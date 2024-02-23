// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2005 - 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check Vector initialization with std::initializer_list objects

#include <deal.II/lac/vector.h>

#include "../tests.h"



int
main()
{
  initlog();

  Vector<double> vd({1.0, 2.0, 3.0});
  for (const auto &x : vd)
    deallog << x << ' ';
  deallog << std::endl;


  Vector<double> vi({1, 2, 3});
  for (const auto &x : vi)
    deallog << x << ' ';
  deallog << std::endl;
}
