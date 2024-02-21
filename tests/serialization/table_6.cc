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


// check serialization for Table<6, int>

#include <deal.II/base/table.h>

#include <boost/serialization/vector.hpp>

#include "serialization.h"


void
test()
{
  unsigned int index1 = 3, index2 = 4, index3 = 2, index4 = 5, index5 = 1,
               index6 = 7;
  TableIndices<6> indices1(index1, index2, index3, index4, index5, index6);
  unsigned int    sum_of_indices =
    index1 + index2 + index3 + index4 + index5 + index6;

  Table<6, int> t1(index1, index2, index3, index4, index5, index6);
  Table<6, int> t2(index1, index2, index3, index4, index5, index6);

  index1 = 2;
  index2 = 5;
  index3 = 4;
  index4 = 1;
  index5 = 5;
  index6 = 8;
  Table<6, int> t3(index1, index2, index3, index4, index5, index6);

  unsigned int counter = 0;
  for (unsigned int i1 = 0; i1 < indices1[0]; ++i1)
    {
      for (unsigned int i2 = 0; i2 < indices1[1]; ++i2)
        {
          for (unsigned int i3 = 0; i3 < indices1[2]; ++i3)
            {
              for (unsigned int i4 = 0; i4 < indices1[3]; ++i4)
                {
                  for (unsigned int i5 = 0; i5 < indices1[4]; ++i5)
                    {
                      for (unsigned int i6 = 0; i6 < indices1[5]; ++i6)
                        {
                          t1[i1][i2][i3][i4][i5][i6] = counter++;
                          t2[i1][i2][i3][i4][i5][i6] = counter + sum_of_indices;
                        }
                    }
                }
            }
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
