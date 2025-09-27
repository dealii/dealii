// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

//
// Verify that MultithreadInfo::set_thread_limit works as advertised
//

#include <deal.II/base/multithread_info.h>

#ifdef DEAL_II_TBB_WITH_ONEAPI
#  include <oneapi/tbb/global_control.h>
#  include <tbb/info.h>
#endif

#include "../tests.h"

int
main()
{
  initlog();

#ifdef DEAL_II_TBB_WITH_ONEAPI
  //
  // The thread concurrency of testsuite executables can be limited by
  // various means. In particular, by design we cannot increase thread
  // concurrency past the environment variable DEAL_II_NUM_THREADS.
  //
  // We thus do the following: We start with the current number of threads
  // - which is a known good value - and iteratively decrease the number by
  // one and check that the tbb thread pool adjusts accordingly:
  //

  for (unsigned int n = MultithreadInfo::n_threads(); n > 0; --n)
    {
      dealii::MultithreadInfo::set_thread_limit(n);

      const auto n_threads = dealii::MultithreadInfo::n_threads();
      const auto n_tbb     = tbb::global_control::active_value(
        tbb::global_control::max_allowed_parallelism);

      if (n != n_threads || n != n_tbb)
        deallog
          << "Problem: Thread limits differ from what has been enforced:\n"
          << "dealii::MultiThreadInfo::set_thread_limit(" << n << ")\n"
          << "dealii::MultiThreadInfo::n_threads() == " << n_threads << "\n"
          << "tbb::global_control::active_value(...) == " << n_tbb << "\n"
          << std::endl;
    }
#endif

  deallog << "OK" << std::endl;
}
