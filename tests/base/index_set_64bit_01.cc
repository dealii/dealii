// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test IndexSet with 64 bit indices

#include <deal.II/base/index_set.h>

#include "../tests.h"


void
print(const IndexSet &s)
{
  deallog << "IndexSet: size " << s.size() << "\n"
          << "n_elements " << s.n_elements() << "\n"
          << "n_intervals " << s.n_intervals() << "\n"
          << "largest_range_starting_index " << s.largest_range_starting_index()
          << "\n"
          << "is_element(1) " << s.is_element(1) << "\n";

  deallog << "contents:" << std::endl;
  s.print(deallog.get_file_stream());
}


void
test()
{
  const unsigned long long base   = 10000000000ull; // 10 billion
  const unsigned long long N      = 5 * base;
  const unsigned long long large  = 2 * base;
  const unsigned long long large2 = 3 * base;
  const unsigned long long large3 = 4 * base;

  {
    deallog << "* complete:" << std::endl;
    IndexSet s = complete_index_set(N);
    Assert(s.size() == N, ExcInternalError());
    Assert(s.n_elements() == N, ExcInternalError());
    print(s);

    IndexSet::size_type x = s.index_within_set(large2);
    deallog << "global: " << large2 << " index_within_set: " << x << std::endl;
    IndexSet::size_type y = s.nth_index_in_set(large2);
    deallog << "local: " << large2 << " nth_index_in_set: " << y << std::endl;
  }
  {
    deallog << "* interval far in:" << std::endl;
    IndexSet s(N);
    s.add_range(large, large2);
    Assert(s.size() == N, ExcInternalError());
    print(s);

    IndexSet::size_type x = s.index_within_set(large + 1);
    deallog << "global: " << large + 1 << " index_within_set: " << x
            << std::endl;
    IndexSet::size_type y = s.nth_index_in_set(3);
    deallog << "local: " << 3 << " nth_index_in_set: " << y << std::endl;
  }
  {
    deallog << "get view:" << std::endl;
    IndexSet s(N);
    s.add_range(large, large2);
    s = s.get_view(large - 1, large2 - 1);
    print(s);
  }
  {
    deallog << "* complete minus interval:" << std::endl;
    IndexSet s = complete_index_set(N);
    IndexSet other(N);
    other.add_range(large, large2);
    s.subtract_set(other);

    Assert(s.size() == N, ExcInternalError());
    print(s);

    IndexSet::size_type x = s.index_within_set(large2 + 1);
    deallog << "global: " << large2 + 1 << " index_within_set: " << x
            << std::endl;
    IndexSet::size_type y = s.nth_index_in_set(x);
    deallog << "local: " << x << " nth_index_in_set: " << y << std::endl;
  }

  {
    deallog << "* add_indices far in:" << std::endl;
    IndexSet s(N);
    IndexSet other(10);
    other.add_index(1);
    other.add_range(4, 7);
    s.add_indices(other, large);
    Assert(s.size() == N, ExcInternalError());
    print(s);
  }

  {
    deallog << "* two large ranges:" << std::endl;
    IndexSet s(N);
    s.add_range(1, large / 2);
    s.add_range(large, large3);
    print(s);

    IndexSet::size_type x = s.index_within_set(large2);
    deallog << "global: " << large2 << " index_within_set: " << x << std::endl;
    IndexSet::size_type y = s.nth_index_in_set(x);
    deallog << "local: " << x << " nth_index_in_set: " << y << std::endl;
  }
}



int
main()
{
  initlog();

  test();
}
