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



// test that DynamicSparsityPattern::iterator yields something
// sensible for empty sparsity patterns

#include <deal.II/lac/dynamic_sparsity_pattern.h>

#include "../tests.h"


void
test()
{
  DynamicSparsityPattern sp(5, 5);
  // put nothing into the sparsity pattern
  sp.compress();

  DynamicSparsityPattern::const_iterator i = sp.begin();
  for (; i != sp.end(); ++i)
    deallog << i->row() << ' ' << i->column() << std::endl;

  deallog << "OK" << std::endl;
}



int
main()
{
  initlog();

  test();
}
