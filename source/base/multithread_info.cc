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

#include <deal.II/base/exceptions.h>
#include <deal.II/base/multithread_info.h>
#include <deal.II/base/types.h>
#include <deal.II/base/utilities.h>

#include <algorithm>
#include <cstdlib> // for std::getenv
#include <mutex>
#include <thread>

#ifdef DEAL_II_WITH_TBB
#  ifdef DEAL_II_TBB_WITH_ONEAPI
#    include <tbb/global_control.h>
#  else
#    include <tbb/task_scheduler_init.h>
#  endif
#endif


#ifdef DEAL_II_WITH_TASKFLOW
#  include <taskflow/taskflow.hpp>
#endif

DEAL_II_NAMESPACE_OPEN


unsigned int
MultithreadInfo::n_cores()
{
  // There is a slight semantic change between our n_cores() call and the
  // std::thread alternative: in case of an error the latter one returns 0
  // in contrast to a 1 that n_cores() used to do. For compatibility, let's
  // translate to our numbering scheme:
  const unsigned int n_cores = std::thread::hardware_concurrency();
  return n_cores == 0 ? 1 : n_cores;
}


void
MultithreadInfo::set_thread_limit(const unsigned int max_threads)
{
  // set the maximal number of threads to the given value as specified
  n_max_threads = max_threads;

  // then also see if something was given in the environment
  {
    if (const char *penv = std::getenv("DEAL_II_NUM_THREADS"))
      {
        unsigned int max_threads_env = numbers::invalid_unsigned_int;
        try
          {
            max_threads_env = Utilities::string_to_int(std::string(penv));
          }
        catch (...)
          {
            AssertThrow(
              false,
              ExcMessage(
                std::string(
                  "When specifying the <DEAL_II_NUM_THREADS> environment "
                  "variable, it needs to be something that can be interpreted "
                  "as an integer. The text you have in the environment "
                  "variable is <") +
                penv + ">"));
          }

        AssertThrow(max_threads_env > 0,
                    ExcMessage(
                      "When specifying the <DEAL_II_NUM_THREADS> environment "
                      "variable, it needs to be a positive number."));

        if (n_max_threads != numbers::invalid_unsigned_int)
          n_max_threads = std::min(n_max_threads, max_threads_env);
        else
          n_max_threads = max_threads_env;
      }
  }

  // If we have not set the number of allowed threads yet, just default to
  // the number of available cores
  if (n_max_threads == numbers::invalid_unsigned_int)
    n_max_threads = n_cores();

#ifdef DEAL_II_WITH_TBB
#  ifdef DEAL_II_TBB_WITH_ONEAPI
  // tbb::global_control is a class that affects the specified behavior of
  // tbb during its lifetime. Thus, in order to set a global thread limit
  // for tbb we have to maintain the object throughout the execution of the
  // program. We do this by maintaining a static std::unique_ptr.
  //
  // A std::unique_ptr is a good choice here because tbb::global_control
  // does not provide a mechanism to override its setting - we can only
  // delete the old and replace it with a new one.
  static std::unique_ptr<tbb::global_control> tbb_global_control;
  tbb_global_control = std::make_unique<tbb::global_control>(
    tbb::global_control::max_allowed_parallelism, n_max_threads);

#  else
  // Initialize the scheduler and destroy the old one before doing so
  static tbb::task_scheduler_init dummy(tbb::task_scheduler_init::deferred);
  if (dummy.is_active())
    dummy.terminate();
  dummy.initialize(n_max_threads);
#  endif
#endif

#ifdef DEAL_II_WITH_TASKFLOW
  executor = std::make_unique<tf::Executor>(n_max_threads);
#endif
}



unsigned int
MultithreadInfo::n_threads()
{
  Assert(n_max_threads != numbers::invalid_unsigned_int, ExcInternalError());
  return n_max_threads;
}



bool
MultithreadInfo::is_running_single_threaded()
{
  return n_threads() == 1;
}



std::size_t
MultithreadInfo::memory_consumption()
{
  // only simple data elements, so use sizeof operator
  return sizeof(MultithreadInfo);
}



void
MultithreadInfo::initialize_multithreading()
{
  static std::once_flag is_initialized;
  std::call_once(is_initialized, []() {
    MultithreadInfo::set_thread_limit(numbers::invalid_unsigned_int);
  });
}



#ifdef DEAL_II_WITH_TASKFLOW
tf::Executor &
MultithreadInfo::get_taskflow_executor()
{
  // This should not trigger in normal user code, because we initialize the
  // Executor in the static DoOnce struct at the end of this file unless you
  // ask for the Executor before this static object gets constructed.
  Assert(
    executor.get() != nullptr,
    ExcMessage(
      "Please initialize multithreading using MultithreadInfo::set_thread_limit() first."));
  return *(executor.get());
}

std::unique_ptr<tf::Executor> MultithreadInfo::executor = nullptr;
#endif

unsigned int MultithreadInfo::n_max_threads = numbers::invalid_unsigned_int;

namespace
{
  // Force the first call to set_thread_limit happen before any tasks in TBB are
  // used. This is necessary as tbb::task_scheduler_init has no effect if TBB
  // got automatically initialized (which happens the first time we use it).
  struct DoOnce
  {
    DoOnce()
    {
      MultithreadInfo::initialize_multithreading();
    }
  } do_once;
} // namespace

DEAL_II_NAMESPACE_CLOSE
