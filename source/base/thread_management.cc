// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2017 by the deal.II authors
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
    static std::atomic<unsigned int> n_existing_threads_counter(1);


    void
    register_thread()
    {
      ++n_existing_threads_counter;
    }



    void
    deregister_thread()
    {
      --n_existing_threads_counter;
      Assert(n_existing_threads_counter >= 1, ExcInternalError());
    }



    [[noreturn]] void
    handle_std_exception(const std::exception &exc) {
      // lock the following context
      // to ensure that we don't
      // print things over each other
      // if we have trouble from
      // multiple threads. release
      // the lock before calling
      // std::abort, though
      static Mutex mutex;
      {
        Mutex::ScopedLock lock(mutex);

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



      [[noreturn]] void handle_unknown_exception()
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
        Mutex::ScopedLock lock(mutex);

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



  unsigned int
  n_existing_threads()
  {
    return internal::n_existing_threads_counter;
  }


  unsigned int
  this_thread_id()
  {
#ifdef SYS_gettid
    const pid_t this_id = syscall(SYS_gettid);
#elif defined(DEAL_II_HAVE_UNISTD_H) && defined(DEAL_II_HAVE_GETPID)
    const pid_t this_id = getpid();
#else
    const unsigned int this_id = 0;
#endif

    return static_cast<unsigned int>(this_id);
  }



#ifdef DEAL_II_USE_MT_POSIX


#  ifndef DEAL_II_USE_MT_POSIX_NO_BARRIERS
  PosixThreadBarrier::PosixThreadBarrier(const unsigned int count,
                                         const char *,
                                         void *)
  {
    pthread_barrier_init(&barrier, nullptr, count);
  }

#  else

  PosixThreadBarrier::PosixThreadBarrier(const unsigned int count,
                                         const char *,
                                         void *)
    : count(count)
  {
    // throw an exception unless we
    // have the special case that a
    // count of 1 is given, since
    // then waiting for a barrier is
    // a no-op, and we don't need the
    // POSIX functionality
    AssertThrow(count == 1,
                ExcMessage("Your local POSIX installation does not support\n"
                           "POSIX barriers. You will not be able to use\n"
                           "this class, but the rest of the threading\n"
                           "functionality is available."));
  }
#  endif



  PosixThreadBarrier::~PosixThreadBarrier()
  {
#  ifndef DEAL_II_USE_MT_POSIX_NO_BARRIERS
    pthread_barrier_destroy(&barrier);
#  else
    // unless the barrier is a no-op,
    // complain again (how did we get
    // here then?)
    if (count != 1)
      std::abort();
#  endif
  }



  int
  PosixThreadBarrier::wait()
  {
#  ifndef DEAL_II_USE_MT_POSIX_NO_BARRIERS
    return pthread_barrier_wait(&barrier);
#  else
    // in the special case, this
    // function is a no-op. otherwise
    // complain about the missing
    // POSIX functions
    if (count == 1)
      return 0;
    else
      {
        std::abort();
        return 1;
      };
#  endif
  }



#endif



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
      };
    return return_values;
  }
} // namespace Threads


DEAL_II_NAMESPACE_CLOSE
