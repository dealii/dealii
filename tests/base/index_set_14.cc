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

#include <stdlib.h>

#include "../tests.h"


void
test()
{
  IndexSet is1(100);
  IndexSet is2(100);

  // randomly add 90 elements to each
  // set, some of which may be
  // repetitions of previous ones
  for (unsigned int i = 0; i < 9 * is1.size() / 10; ++i)
    {
      is1.add_index(Testing::rand() % is1.size());
      is2.add_index(Testing::rand() % is2.size());
    }

  IndexSet is3 = is1 & is2;

  deallog << "Set sizes: " << is1.n_elements() << ' ' << is2.n_elements() << ' '
          << is3.n_elements() << std::endl;

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
