// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2004 - 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------



// test that i++ and i+n operations work with SparsityPattern::iterator.

#include <deal.II/lac/sparsity_pattern.h>

#include "../tests.h"

void
test()
{
  SparsityPattern sp(5, 5, 3);
  for (unsigned int i = 0; i < 5; ++i)
    for (unsigned int j = 0; j < 5; ++j)
      if ((i + 2 * j + 1) % 3 == 0)
        sp.add(i, j);
  sp.compress();

  SparsityPattern::const_iterator i = sp.begin();

  ++i;
  i++;
  i + 2;

  deallog << "OK" << std::endl;
}



int
main()
{
  initlog();

  test();
}
