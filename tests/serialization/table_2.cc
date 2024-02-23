// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// check serialization for Table<2, int>

#include <deal.II/base/table.h>

#include <boost/serialization/vector.hpp>

#include "serialization.h"


void
test()
{
  unsigned int    index1 = 3, index2 = 4;
  TableIndices<2> indices1(index1, index2);
  unsigned int    sum_of_indices = index1 + index2;

  Table<2, int> t1(index1, index2);
  Table<2, int> t2(index1, index2);

  index1 = 2;
  index2 = 5;
  Table<2, int> t3(index1, index2);

  unsigned int counter = 0;
  for (unsigned int i1 = 0; i1 < indices1[0]; ++i1)
    {
      for (unsigned int i2 = 0; i2 < indices1[1]; ++i2)
        {
          t1[i1][i2] = counter++;
          t2[i1][i2] = counter + sum_of_indices;
        }
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
