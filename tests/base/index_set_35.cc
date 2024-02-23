// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Ensure that we can compare IndexSet objects against the empty
// object. (Comparisons between IndexSet objects of different sizes
// are otherwise not allowed.)

#include <deal.II/base/index_set.h>

#include <stdlib.h>

#include "../tests.h"


void
test()
{
  IndexSet is1(100);
  IndexSet is2;

  Assert((is1 == is2) == false, ExcInternalError());
  Assert((is2 != is1) == true, ExcInternalError());
  Assert((is1 == is2) == false, ExcInternalError());
  Assert((is2 != is1) == true, ExcInternalError());

  deallog << "OK" << std::endl;
}



int
main()
{
  initlog();

  test();
}
