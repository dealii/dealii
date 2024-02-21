// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test IndexSet iterators

#include <deal.II/base/index_set.h>

#include "../tests.h"

void
test(IndexSet &index_set)
{
  index_set.print(deallog);

  Assert((int)index_set.n_intervals() ==
           index_set.end_intervals() - index_set.begin_intervals(),
         ExcInternalError());

  IndexSet::IntervalIterator endit = index_set.end_intervals();
  Assert(!endit->is_valid(), ExcInternalError());

  // print intervals
  for (IndexSet::IntervalIterator it = index_set.begin_intervals(); it != endit;
       ++it)
    {
      if (it->n_elements() == 1)
        deallog << *(it->begin()) << ' ';
      else
        deallog << '[' << *(it->begin()) << ',' << it->last() << "] ";
    }
  deallog << std::endl;

  // print entries
  {
    for (IndexSet::ElementIterator it = index_set.begin();
         it != index_set.end();
         ++it)
      deallog << *it << ' ';
    deallog << std::endl;
  }

  // check comparison, distance, and n_elements:
  unsigned int c = 0;
  for (IndexSet::IntervalIterator it = index_set.begin_intervals(); it != endit;
       ++it, ++c)
    {
      IndexSet::IntervalIterator it2 = it;
      Assert(it == it2, ExcInternalError());
      Assert(it - it2 == 0, ExcInternalError());
      Assert(endit != it, ExcInternalError());
      Assert(it != endit, ExcInternalError());

      IndexSet::IntervalIterator it3 = it2++;
      Assert(it == it3, ExcInternalError());
      Assert(it2 - it3 == 1, ExcInternalError());
      Assert(it3 < it2, ExcInternalError());
      Assert(++it3 == it2, ExcInternalError());

      Assert(it->is_valid(), ExcInternalError());

      Assert(it < endit, ExcInternalError());
      Assert(!(it < it), ExcInternalError());

      Assert((it - index_set.begin_intervals()) == (int)c, ExcInternalError());
      Assert((index_set.begin_intervals() - it) == -(int)c, ExcInternalError());

      deallog << c << ": n_el: " << it->n_elements() << std::endl;
    }

  // ElementIterator checks
  {
    IndexSet::ElementIterator it = index_set.begin();

    unsigned int c = 0;
    for (; it != index_set.end(); ++it, ++c)
      {
        Assert(it < index_set.end(), ExcInternalError());

        IndexSet::ElementIterator it2 = it;
        Assert(it == it2, ExcInternalError());
        Assert(!(it < it2), ExcInternalError());
        IndexSet::ElementIterator it3 = it2++;
        Assert(it == it3, ExcInternalError());
        Assert(it < it2, ExcInternalError());
        Assert(it != it2, ExcInternalError());

        Assert((it - index_set.begin()) == c, ExcInternalError());
        Assert((index_set.begin() - it) == -(int)c, ExcInternalError());
      }
  }


  // pre vs post increment
  {
    IndexSet::ElementIterator it  = index_set.begin();
    IndexSet::ElementIterator it2 = index_set.begin();

    for (; it != index_set.end(); ++it, it2++)
      Assert(it == it2, ExcInternalError());
  }

  // pre vs post increment, part 2
  {
    IndexSet::IntervalIterator it  = index_set.begin_intervals();
    IndexSet::IntervalIterator it2 = index_set.begin_intervals();

    for (; it != index_set.end_intervals(); ++it, it2++)
      Assert(it == it2, ExcInternalError());
  }

  // Assignment
  {
    IndexSet::IntervalIterator it = index_set.begin_intervals();
    Assert(it->is_valid() || index_set.n_elements() == 0, ExcInternalError());
    it = index_set.end_intervals();
    Assert(!it->is_valid(), ExcInternalError());
  }

  // Construct empty iterator
  {
    IndexSet::IntervalIterator it;
    Assert(!it->is_valid(), ExcInternalError());
    it = index_set.begin_intervals();
    Assert(it->is_valid() || index_set.n_elements() == 0, ExcInternalError());
  }
}

void
test()
{
  IndexSet index_set(20);
  index_set.add_range(2, 4);
  index_set.add_range(12, 18);
  index_set.add_index(6);
  index_set.add_index(8);
  index_set.add_index(9);
  index_set.add_index(16);

  test(index_set);

  IndexSet empty(42);
  test(empty);
}



int
main()
{
  initlog();

  test();
}
