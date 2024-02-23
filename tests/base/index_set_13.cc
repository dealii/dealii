// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2009 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test IndexSet::operator &

#include <deal.II/base/index_set.h>

#include "../tests.h"


void
test()
{
  IndexSet is1(20);
  IndexSet is2(20);

  is1.add_index(2);
  is1.add_index(3);
  is1.add_index(4);
  is1.add_index(6);
  is1.add_index(7);

  is2.add_range(4, 9);

  IndexSet is3 = is1 & is2;

  for (unsigned int i = 0; i < is3.size(); ++i)
    {
      deallog << i << ' ' << (is3.is_element(i) ? "true" : "false")
              << std::endl;

      AssertThrow((is1.is_element(i) && is2.is_element(i)) == is3.is_element(i),
                  ExcInternalError());
    }

  // some sanity tests
  AssertThrow((is1 & is2) == (is2 & is1), ExcInternalError());
  AssertThrow((is1 & is3) == (is2 & is3), ExcInternalError());
  AssertThrow((is1 & is3) == is3, ExcInternalError());
  AssertThrow((is3 & is1) == is3, ExcInternalError());
  AssertThrow((is3 & is3) == is3, ExcInternalError());

  deallog << "OK" << std::endl;
}



int
main()
{
  initlog();

  test();
}
