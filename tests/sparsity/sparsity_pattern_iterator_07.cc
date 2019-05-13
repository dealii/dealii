// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2018 by the deal.II authors
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



// test that SparsityPattern::iterator yields something
// sensible for empty sparsity patterns

#include <deal.II/lac/sparsity_pattern.h>

#include "../tests.h"


void
test()
{
  // use a non-square sparsity pattern to make sure the object does
  // not automatically store the diagonals
  SparsityPattern sp(5, 6, 3);
  // put nothing into the sparsity pattern
  sp.compress();

  SparsityPattern::const_iterator i = sp.begin();
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
