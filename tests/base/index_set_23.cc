// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// we managed to get a function as simple as IndexSet::operator== wrong -- who
// knew?

#include <deal.II/base/index_set.h>

#include <stdlib.h>

#include "../tests.h"


void
test()
{
  IndexSet is1(100), is2(100);

  is1.add_range(0, 10);
  is2.add_range(0, 20);

  Assert((is1 == is2) == false, ExcInternalError());
  Assert((is1 != is2) == true, ExcInternalError());

  deallog << "OK" << std::endl;
}



int
main()
{
  initlog();

  test();
}
