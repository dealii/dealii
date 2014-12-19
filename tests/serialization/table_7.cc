// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


// check serialization for Table<7, int>

#include "serialization.h"
#include <deal.II/base/table.h>
#include <boost/serialization/vector.hpp>

void test ()
{
  unsigned int index1 = 3, index2 = 4, index3 = 2, index4 = 5, index5 = 1,
               index6 = 7, index7 = 6;
  TableIndices<7> indices1(index1, index2, index3, index4, index5, index6, index7);
  unsigned int sum_of_indices = index1 + index2 + index3 + index4 + index5
                                + index6 + index7;

  Table<7, int> t1(index1, index2, index3, index4, index5, index6, index7);
  Table<7, int> t2(index1, index2, index3, index4, index5, index6, index7);

  index1 = 2;
  index2 = 5;
  index3 = 4;
  index4 = 1;
  index5 = 5;
  index6 = 8;
  index7 = 6;
  Table<7, int> t3(index1, index2, index3, index4, index5, index6, index7);

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
                          for (unsigned int i7 = 0; i7 < indices1[6]; ++i7)
                            {
                              t1[i1][i2][i3][i4][i5][i6][i7] = counter ++;
                              t2[i1][i2][i3][i4][i5][i6][i7] = counter + sum_of_indices;
                            }
                        }
                    }
                }
            }
        }
    }

  verify (t1, t2);

  verify (t1, t3);
}


int main ()
{
  std::ofstream logfile("output");
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test ();

  deallog << "OK" << std::endl;
}
