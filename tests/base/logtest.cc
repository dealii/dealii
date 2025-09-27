// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 1998 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// some tests for the logstream class, writing output, pushing and
// popping levels, etc


#include <limits>

#include "../tests.h"


int
main()
{
  initlog();

  deallog << "Test" << std::endl;
  deallog.push("l1");
  deallog << "Test1" << std::endl;
  deallog.push("l2");
  deallog << "Test2"
          << "Test3" << std::endl;
  deallog.push("l3");
  deallog << "Test4";
  deallog.pop();
  deallog << "Test5" << std::endl;
  deallog.pop();
  deallog << "Test6" << std::endl;
  deallog.pop();
  deallog << "Test7" << std::endl;
}
