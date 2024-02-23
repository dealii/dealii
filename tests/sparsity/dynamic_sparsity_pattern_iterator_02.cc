// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2004 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// test DynamicSparsityPattern::iterator with sparsity patterns that
// have empty rows

#include <deal.II/lac/dynamic_sparsity_pattern.h>

#include "../tests.h"


void
test()
{
  DynamicSparsityPattern sp(5, 5);
  sp.add(0, 0);
  sp.add(3, 3);
  sp.compress();

  DynamicSparsityPattern::const_iterator i = sp.begin();
  for (; i != sp.end(); ++i)
    deallog << i->row() << ' ' << i->column() << std::endl;

  deallog << "OK" << std::endl;

  i = sp.begin(1);
  deallog << i->row() << ' ' << i->column() << std::endl;
  deallog << "OK" << std::endl;
}



int
main()
{
  initlog();

  test();
}
