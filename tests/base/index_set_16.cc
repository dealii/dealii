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


// test IndexSet::subtract_set further

#include <deal.II/base/index_set.h>

#include <stdlib.h>

#include "../tests.h"


void
test()
{
  IndexSet is1(100);
  IndexSet is2(100);

  is1.add_range(0, 10);
  is1.add_range(20, 100);

  is2.add_range(0, 50);
  //  is2.add_range(10,15);

  IndexSet is3 = is1;
  is3.subtract_set(is2);

  is1.print(deallog);
  is2.print(deallog);
  is3.print(deallog);

  for (unsigned int i = 0; i < is3.size(); ++i)
    {
      AssertThrow((is1.is_element(i) && !is2.is_element(i)) ==
                    is3.is_element(i),
                  ExcInternalError());
    }

  deallog << is3.index_within_set(51) << std::endl;
  AssertThrow(is3.index_within_set(51) == 1, ExcInternalError());

  deallog << "OK" << std::endl;
}



int
main()
{
  initlog();

  test();
}
