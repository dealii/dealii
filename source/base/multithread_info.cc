// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2020 by the deal.II authors
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

#include <deal.II/base/multithread_info.h>
#include <deal.II/base/utilities.h>

#include <algorithm>
#include <cstdlib> // for std::getenv
#include <thread>

#ifdef DEAL_II_WITH_TBB
#  define TBB_SUPPRESS_DEPRECATED_MESSAGES 1
#  include <tbb/task_scheduler_init.h>
#  undef TBB_SUPPRESS_DEPRECATED_MESSAGES
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
  // Initialize the scheduler and destroy the old one before doing so
  static tbb::task_scheduler_init dummy(tbb::task_scheduler_init::deferred);
  if (dummy.is_active())
    dummy.terminate();
  dummy.initialize(n_max_threads);
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
  static bool done = false;
  if (done)
    return;

  MultithreadInfo::set_thread_limit(numbers::invalid_unsigned_int);
  done = true;
}


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
