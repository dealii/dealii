// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2019 by the deal.II authors
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
