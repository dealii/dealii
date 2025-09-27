// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// check serialization for Table<1. int>

#include <deal.II/base/table.h>

#include <boost/serialization/vector.hpp>

#include "serialization.h"

void
test()
{
  unsigned int  index1 = 3;
  Table<1, int> t1(index1);

  Table<1, int> t2(index1);

  unsigned int  index3 = 2;
  Table<1, int> t3(index3);

  for (unsigned int i = 0; i < index1; ++i)
    {
      t1[i] = i + 1;
      t2[i] = i + 1 + index1;
    }
  verify(t1, t2);

  verify(t1, t3);
}


int
main()
{
  initlog();
  deallog << std::setprecision(3);

  test();

  deallog << "OK" << std::endl;
}
