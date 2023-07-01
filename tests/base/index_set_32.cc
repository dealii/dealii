// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2018 by the deal.II authors
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


// Test how IndexSet::add_indices deals with duplicate indices to be added

#include <deal.II/base/index_set.h>

#include "../tests.h"


void
test()
{
  IndexSet index_set(10);

  {
    const unsigned int array[] = {2, 2, 2, 3, 4, 6, 8, 8};
    index_set.add_indices((const unsigned int *)array,
                          array + sizeof(array) / sizeof(array[0]));
  }

  Assert(index_set.is_contiguous() == false, ExcInternalError());

  for (unsigned int i = 0; i < index_set.size(); ++i)
    deallog << i << ' ' << (index_set.is_element(i) ? "true" : "false")
            << std::endl;
}



int
main()
{
  initlog();

  test();
}
