// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2012 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// tests for the BlockMask class
//
// here: test that creating a mask from a vector<bool> works


#include <deal.II/fe/block_mask.h>

#include "../tests.h"



void
test()
{
  std::vector<bool> v(12);
  for (unsigned int i = 0; i < v.size(); ++i)
    v[i] = (i % 3 == 0);

  BlockMask m(v);

  // verify equality
  for (unsigned int i = 0; i < v.size(); ++i)
    Assert(m[i] == v[i], ExcInternalError());

  // this needs to throw an exception
  try
    {
      m[v.size()];
    }
  catch (ExceptionBase &e)
    {
      deallog << e.get_exc_name() << std::endl;
    }

  deallog << "OK" << std::endl;
}


int
main()
{
  deal_II_exceptions::disable_abort_on_exception();

  initlog();
  deallog << std::setprecision(4);

  test();
}
