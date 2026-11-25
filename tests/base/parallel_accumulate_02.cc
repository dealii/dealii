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


// test parallel::accumulate_from_subranges with integral ranges that do not
// start at zero, random access and bidirectional iterator ranges, empty
// ranges, floating point results, and calls from within tasks

#include <deal.II/base/parallel.h>
#include <deal.II/base/thread_management.h>

#include <list>
#include <numeric>
#include <vector>

#include "../tests.h"


void
test(const bool print)
{
  // integral range not starting at zero
  {
    const long long s = parallel::accumulate_from_subranges<long long>(
      [](const std::size_t begin, const std::size_t end) {
        long long s = 0;
        for (std::size_t i = begin; i < end; ++i)
          s += i;
        return s;
      },
      std::size_t(17),
      std::size_t(100017),
      10);
    AssertThrow(s == 100017LL * 100016 / 2 - 17 * 16 / 2, ExcInternalError());
    if (print)
      deallog << "integral: " << s << std::endl;
  }

  // random access iterators and floating point accumulation
  {
    std::vector<double> x(100000);
    for (unsigned int i = 0; i < x.size(); ++i)
      x[i] = 1. / (1. + i);
    const double reference = std::accumulate(x.begin(), x.end(), 0.);

    const double s = parallel::accumulate_from_subranges<double>(
      [](const std::vector<double>::const_iterator begin,
         const std::vector<double>::const_iterator end) {
        return std::accumulate(begin, end, 0.);
      },
      x.cbegin(),
      x.cend(),
      100);
    AssertThrow(std::abs(s - reference) < 1e-12 * reference,
                ExcInternalError());
    if (print)
      deallog << "random access: OK" << std::endl;
  }

  // bidirectional iterators
  {
    const std::list<int> l(5000, 2);
    const int            s = parallel::accumulate_from_subranges<int>(
      [](const std::list<int>::const_iterator begin,
         const std::list<int>::const_iterator end) {
        return std::accumulate(begin, end, 0);
      },
      l.cbegin(),
      l.cend(),
      7);
    AssertThrow(s == 10000, ExcInternalError());
    if (print)
      deallog << "bidirectional: " << s << std::endl;
  }

  // range smaller than the grain size, and an empty range
  {
    const auto length = [](const int begin, const int end) {
      return end - begin;
    };
    AssertThrow(parallel::accumulate_from_subranges<int>(length, 0, 5, 100) ==
                  5,
                ExcInternalError());
    AssertThrow(parallel::accumulate_from_subranges<int>(length, 3, 3, 1) == 0,
                ExcInternalError());
    if (print)
      deallog << "small and empty ranges: OK" << std::endl;
  }
}


int
main()
{
  initlog();

  test(true);

  // Also call the function concurrently from within several tasks:
  Threads::TaskGroup<void> tasks;
  for (unsigned int i = 0; i < 4; ++i)
    tasks += Threads::new_task([]() { test(false); });
  tasks.join_all();

  deallog << "OK" << std::endl;
}
