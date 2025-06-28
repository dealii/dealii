// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2000 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/mutex.h>
#include <deal.II/base/thread_management.h>

#include <atomic>
#include <cstdlib>
#include <iostream>

#ifdef DEAL_II_HAVE_UNISTD_H
#  include <unistd.h>
#endif


DEAL_II_NAMESPACE_OPEN


namespace Threads
{
  namespace internal
  {
    [[noreturn]] void
    handle_std_exception(const std::exception &exc)
    {
      // lock the following context
      // to ensure that we don't
      // print things over each other
      // if we have trouble from
      // multiple threads. release
      // the lock before calling
      // std::abort, though
      static Mutex mutex;
      {
        std::lock_guard<std::mutex> lock(mutex);

        std::cerr
          << std::endl
          << std::endl
          << "---------------------------------------------------------"
          << std::endl
          << "In one of the sub-threads of this program, an exception\n"
          << "was thrown and not caught. Since exceptions do not\n"
          << "propagate to the main thread, the library has caught it.\n"
          << "The information carried by this exception is given below.\n"
          << std::endl
          << "---------------------------------------------------------"
          << std::endl;
        std::cerr << "Exception message: " << std::endl
                  << "  " << exc.what() << std::endl
                  << "Exception type: " << std::endl
                  << "  " << typeid(exc).name() << std::endl;
        std::cerr << "Aborting!" << std::endl
                  << "---------------------------------------------------------"
                  << std::endl;
      }

      std::abort();
    }



    [[noreturn]] void
    handle_unknown_exception()
    {
      // lock the following context
      // to ensure that we don't
      // print things over each other
      // if we have trouble from
      // multiple threads. release
      // the lock before calling
      // std::abort, though
      static Mutex mutex;
      {
        std::lock_guard<std::mutex> lock(mutex);

        std::cerr
          << std::endl
          << std::endl
          << "---------------------------------------------------------"
          << std::endl
          << "In one of the sub-threads of this program, an exception\n"
          << "was thrown and not caught. Since exceptions do not\n"
          << "propagate to the main thread, the library has caught it.\n"
          << std::endl
          << "---------------------------------------------------------"
          << std::endl;
        std::cerr << "Type of exception is unknown, but not std::exception.\n"
                  << "No additional information is available.\n"
                  << "---------------------------------------------------------"
                  << std::endl;
      }
      std::abort();
    }
  } // namespace internal



  std::vector<std::pair<unsigned int, unsigned int>>
  split_interval(const unsigned int begin,
                 const unsigned int end,
                 const unsigned int n_intervals)
  {
    Assert(end >= begin, ExcInternalError());

    const unsigned int n_elements              = end - begin;
    const unsigned int n_elements_per_interval = n_elements / n_intervals;
    const unsigned int residual                = n_elements % n_intervals;

    std::vector<std::pair<unsigned int, unsigned int>> return_values(
      n_intervals);

    return_values[0].first = begin;
    for (unsigned int i = 0; i < n_intervals; ++i)
      {
        if (i != n_intervals - 1)
          {
            return_values[i].second =
              (return_values[i].first + n_elements_per_interval);
            // distribute residual in
            // division equally among
            // the first few
            // subintervals
            if (i < residual)
              ++return_values[i].second;
            return_values[i + 1].first = return_values[i].second;
          }
        else
          return_values[i].second = end;
      }
    return return_values;
  }
} // namespace Threads


DEAL_II_NAMESPACE_CLOSE
