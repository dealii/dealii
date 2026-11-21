// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


// With Taskflow, parallel::apply_to_subranges() should only split ranges of
// random access iterators. With other iterators, even ones that have an
// operator+ like those of SparseMatrix, every subrange would have to step over
// all elements before it. Instead, the function should be called once on the
// whole range.

#include <deal.II/base/parallel.h>
#include <deal.II/base/std_cxx20/iota_view.h>

#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparsity_pattern.h>

#include <atomic>
#include <list>
#include <vector>

#include "../tests.h"


template <typename Iterator>
void
test(const Iterator &begin, const Iterator &end)
{
  std::atomic<unsigned int> n_calls(0);
  std::atomic<bool>         whole_range(false);
  parallel::apply_to_subranges(
    begin,
    end,
    [&](const Iterator &b, const Iterator &e) {
      ++n_calls;
      whole_range = (b == begin && e == end);
    },
    /* grainsize = */ 10);

  deallog << "number of calls = " << n_calls.load() << std::endl;
  AssertThrow(whole_range, ExcInternalError());
}


// Ranges of random access iterators should still be split, with every element
// visited exactly once. In C++20, the iterators of std::ranges::iota_view
// report input_iterator_tag but satisfy std::random_access_iterator.
void
test_split()
{
  const auto range =
    std_cxx20::ranges::iota_view<unsigned int, unsigned int>(17, 10017);
  using Iterator = decltype(range.begin());

  std::vector<std::atomic<unsigned int>> n_visits(10000);
  std::atomic<unsigned int>              n_calls(0);
  parallel::apply_to_subranges(
    range.begin(),
    range.end(),
    [&](const Iterator &b, const Iterator &e) {
      ++n_calls;
      for (Iterator i = b; i != e; ++i)
        ++n_visits[*i - 17];
    },
    /* grainsize = */ 10);

  for (const auto &n : n_visits)
    AssertThrow(n == 1, ExcInternalError());
  if (MultithreadInfo::n_threads() > 1)
    AssertThrow(n_calls > 1, ExcInternalError());

  deallog << "OK" << std::endl;
}


int
main()
{
  initlog();

  std::list<unsigned int> list(10000);
  test(list.begin(), list.end());

  SparsityPattern sparsity(10000, 10000, 1);
  sparsity.compress();
  SparseMatrix<double> matrix(sparsity);
  test(matrix.begin(), matrix.end());

  test_split();
}
