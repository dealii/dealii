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


// test parallel::ParallelForInteger::apply_parallel, both when called
// directly and from within tasks

#include <deal.II/base/parallel.h>
#include <deal.II/base/thread_management.h>

#include <vector>

#include "../tests.h"


struct Square : public parallel::ParallelForInteger
{
  Square(std::vector<double> &values)
    : values(values)
  {}

  void
  apply_to_subrange(const std::size_t begin,
                    const std::size_t end) const override
  {
    for (std::size_t i = begin; i < end; ++i)
      values[i] = static_cast<double>(i) * i;
  }

  std::vector<double> &values;
};



void
test(const std::size_t begin,
     const std::size_t end,
     const std::size_t grain_size,
     const bool        print)
{
  std::vector<double> values(end + 3, -1.);

  const Square square(values);
  square.apply_parallel(begin, end, grain_size);

  for (std::size_t i = 0; i < values.size(); ++i)
    AssertThrow(values[i] ==
                  ((i >= begin && i < end) ? static_cast<double>(i) * i : -1.),
                ExcInternalError());

  if (print)
    deallog << "[" << begin << "," << end << ") with grain size " << grain_size
            << ": OK" << std::endl;
}



int
main()
{
  initlog();

  test(0, 100000, 1, true);
  test(13, 100000, 1000, true);
  test(0, 10, 100, true);
  test(5, 5, 1, true);

  // Also call the function concurrently from within several tasks:
  Threads::TaskGroup<void> tasks;
  for (unsigned int i = 0; i < 4; ++i)
    tasks += Threads::new_task([]() { test(7, 100000, 10, false); });
  tasks.join_all();

  deallog << "OK" << std::endl;
}
