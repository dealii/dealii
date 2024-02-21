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


// test IndexSet::is_contiguous and compress()

#include <deal.II/base/index_set.h>

#include "../tests.h"


void
test()
{
  IndexSet index_set(10);
  index_set.add_index(2);
  index_set.add_index(3);
  index_set.add_index(4);

  index_set.add_index(6);
  index_set.add_index(7);

  deallog << (index_set.is_contiguous() ? "true" : "false") << std::endl;
  AssertThrow(index_set.is_contiguous() == false, ExcInternalError());

  for (unsigned int i = 0; i < index_set.size(); ++i)
    deallog << i << ' ' << (index_set.is_element(i) ? "true" : "false")
            << std::endl;
}



int
main()
{
  initlog();

  test();
}
