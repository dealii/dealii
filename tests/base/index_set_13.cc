// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2017 by the deal.II authors
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


// test IndexSet::operator &

#include <deal.II/base/index_set.h>

#include "../tests.h"


void
test()
{
  IndexSet is1(20);
  IndexSet is2(20);

  is1.add_index(2);
  is1.add_index(3);
  is1.add_index(4);
  is1.add_index(6);
  is1.add_index(7);

  is2.add_range(4, 9);

  IndexSet is3 = is1 & is2;

  for (unsigned int i = 0; i < is3.size(); ++i)
    {
      deallog << i << ' ' << (is3.is_element(i) ? "true" : "false")
              << std::endl;

      AssertThrow((is1.is_element(i) && is2.is_element(i)) == is3.is_element(i),
                  ExcInternalError());
    }

  // some sanity tests
  AssertThrow((is1 & is2) == (is2 & is1), ExcInternalError());
  AssertThrow((is1 & is3) == (is2 & is3), ExcInternalError());
  AssertThrow((is1 & is3) == is3, ExcInternalError());
  AssertThrow((is3 & is1) == is3, ExcInternalError());
  AssertThrow((is3 & is3) == is3, ExcInternalError());

  deallog << "OK" << std::endl;
}



int
main()
{
  initlog();

  test();
}
