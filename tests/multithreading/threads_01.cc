// ---------------------------------------------------------------------
//
// Copyright (C) 2013 - 2021 by the deal.II authors
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
