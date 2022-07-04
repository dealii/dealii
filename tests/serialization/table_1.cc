// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2021 by the deal.II authors
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
