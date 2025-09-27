// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2012 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// tests for the BlockMask class
//
// here: test that creating a mask that's empty is always true


#include <deal.II/fe/block_mask.h>

#include "../tests.h"



void
test()
{
  BlockMask m;
  AssertThrow(m[0] == true, ExcInternalError());
  AssertThrow(m[42] == true, ExcInternalError());
  AssertThrow(m[1000000000] == true, ExcInternalError());

  deallog << "OK" << std::endl;
}


int
main()
{
  initlog();
  deallog << std::setprecision(4);

  test();
}
