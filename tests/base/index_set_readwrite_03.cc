// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2022 by the deal.II authors
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


// document a bug in IndexSet::block_read() and block_write(),
// because largest_range was not serialized/reset.

#include <deal.II/base/index_set.h>

#include <stdlib.h>

#include "../tests.h"


void
test()
{
  IndexSet is1(250);
  is1.add_range(125, 200);
  is1.add_range(0, 75);

  {
    std::ofstream out("a.idxset");
    is1.block_write(out);
  }

  IndexSet      is2;
  std::ifstream in("a.idxset");
  is2.block_read(in);

  Assert(is1 == is2, ExcInternalError());

  deallog << is1.is_element(4) << ' ' << is2.is_element(4) << std::endl;

  deallog << "OK" << std::endl;

  std::remove("a.idxset");
}



int
main()
{
  initlog();

  test();
}
