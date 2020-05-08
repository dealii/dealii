// ---------------------------------------------------------------------
//
// Copyright (C) 2012 - 2020 by the deal.II authors
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



// tests for the BlockMask class
//
// here: test that creating a mask with constant elements using the direct
// constructor


#include <deal.II/fe/block_mask.h>

#include "../tests.h"



void
test()
{
  std::vector<bool> v(12, false);
  BlockMask         m(12, false);

  // verify equality
  for (unsigned int i = 0; i < v.size(); ++i)
    AssertThrow(m[i] == v[i], ExcInternalError());

  deallog << "OK" << std::endl;
}


int
main()
{
  initlog();
  deallog << std::setprecision(4);

  test();
}
