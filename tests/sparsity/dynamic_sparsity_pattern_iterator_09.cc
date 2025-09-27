// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2004 - 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// test that DynamicSparsityPattern::iterator can be used
// inside STL containers (i.e. check default constructor)

#include <deal.II/lac/dynamic_sparsity_pattern.h>

#include "../tests.h"


void
test()
{
  DynamicSparsityPattern sp;
  // put nothing into the sparsity pattern
  sp.compress();

  std::vector<std::pair<DynamicSparsityPattern::const_iterator,
                        DynamicSparsityPattern::const_iterator>>
    iterators(1);
  iterators[0].first  = sp.begin();
  iterators[0].second = sp.end();

  deallog << "OK" << std::endl;
}



int
main()
{
  initlog();

  test();
}
