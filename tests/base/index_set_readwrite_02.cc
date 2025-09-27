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


// test IndexSet::block_read() and block_write()

#include <deal.II/base/index_set.h>

#include <stdlib.h>

#include "../tests.h"


void
test()
{
  IndexSet is1(100);

  is1.add_range(0, 10);
  is1.add_range(20, 100);

  {
    std::ofstream out("a.idxset");
    is1.block_write(out);
  }

  IndexSet      is2;
  std::ifstream in("a.idxset");
  is2.block_read(in);

  Assert(is1 == is2, ExcInternalError());


  IndexSet is3(11);
  is3.add_range(3, 5);
  std::ifstream in2("a.idxset");
  is3.block_read(in2);

  Assert(is1 == is3, ExcInternalError());

  deallog << "OK" << std::endl;

  std::remove("a.idxset");
}



int
main()
{
  initlog();

  test();
}
