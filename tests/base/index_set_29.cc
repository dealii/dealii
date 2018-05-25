// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2017 by the deal.II authors
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


// test IndexSet::index_within_set () for global indices which are not
// part of a non-contiguous index set.
// This test has exactly the same output as index_set_12

#include <deal.II/base/index_set.h>

#include "../tests.h"


void
test()
{
  IndexSet index_set(20);

  index_set.add_index(2);
  index_set.add_index(3);
  index_set.add_index(4);

  index_set.add_index(6);
  index_set.add_index(7);

  index_set.compress();

  index_set.add_index(9);

  for (IndexSet::size_type i = 0; i < index_set.n_elements(); ++i)
    {
      deallog << index_set.nth_index_in_set(i) << std::endl;

      AssertThrow(
        index_set.index_within_set(index_set.nth_index_in_set(i) == i),
        ExcInternalError());
    }
  deallog << "OK" << std::endl;

  for (IndexSet::size_type i = 0; i < index_set.size(); ++i)
    {
      const IndexSet::size_type i_out  = index_set.index_within_set(i);
      const bool                within = (i_out != numbers::invalid_dof_index);
      AssertThrow(within == index_set.is_element(i), ExcInternalError());

      if (within)
        deallog << i << ' ' << i_out << std::endl;
    }
}



int
main()
{
  initlog();

  test();
}
