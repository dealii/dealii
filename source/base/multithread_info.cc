// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#include <deal.II/base/multithread_info.h>
#include <deal.II/base/utilities.h>

#ifdef HAVE_UNISTD_H
#  include <unistd.h>
#endif

#if defined(__MACH__) && defined(__APPLE__)
#  include <sys/types.h>
#  include <sys/sysctl.h>
#endif

#if defined(__FreeBSD__)
#  include <stdlib.h>
#endif

#ifdef DEAL_II_WITH_THREADS
#  include <deal.II/base/thread_management.h>
#  include <tbb/task_scheduler_init.h>
#endif

DEAL_II_NAMESPACE_OPEN

#ifdef DEAL_II_WITH_THREADS

/* Detecting how many processors a given machine has is something that
   varies greatly between operating systems. For a few operating
   systems, we have figured out how to do that below, but some others
   are still missing. If you find a way to do this on your favorite
   system, please let us know.
 */


#  if defined(__linux__) ||  defined(__sun__) || defined(__osf__) || defined(_AIX)

unsigned int MultithreadInfo::get_n_cpus()
{
  return sysconf(_SC_NPROCESSORS_ONLN);
}

#  elif defined(__MACH__) && defined(__APPLE__)
// This is only tested on a dual G5 2.5GHz running MacOSX 10.3.6
// and on an Intel Mac Book Pro.
// If it doesn't work please contact the mailinglist.
unsigned int MultithreadInfo::get_n_cpus()
{
  int mib[2];
  int n_cpus;
  size_t len;

  mib[0] = CTL_HW;
  mib[1] = HW_NCPU;
  len = sizeof(n_cpus);
  sysctl(mib, 2, &n_cpus, &len, NULL, 0);

  return n_cpus;
}

#  else

// If you get n_cpus=1 although you are on a multi-processor machine,
// then this may have two reasons: either because the system macros,
// e.g.__linux__, __sgi__, etc. weren't defined by the compiler or the
// detection of processors is really not implemented for your specific
// system. In the first case you can add e.g. -D__sgi__ to your
// compiling flags, in the latter case you need to implement the
// get_n_cpus() function for your system.
//
// In both cases, this #else case is compiled, a fact that you can
// easily verify by uncommenting the following #error directive,
// recompiling and getting a compilation error right at that line.
// After definition of the system macro or the implementation of the
// new detection this #error message during compilation shouldn't
// occur any more.
//
// Please send all new implementations of detection of processors to
// the deal.II mailing list, such that it can be included into the
// next deal.II release.

//#error Detection of Processors not supported on this OS. Setting n_cpus=1 by default.

unsigned int MultithreadInfo::get_n_cpus()
{
  return 1;
}

#  endif

void MultithreadInfo::set_thread_limit(const unsigned int max_threads)
{
  Assert(n_max_threads==numbers::invalid_unsigned_int,
         ExcMessage("Calling set_thread_limit() more than once is not supported!"));

  unsigned int max_threads_env = numbers::invalid_unsigned_int;
  char *penv;
  penv = getenv ("DEAL_II_NUM_THREADS");

  if (penv!=NULL)
    max_threads_env = Utilities::string_to_int(std::string(penv));

  n_max_threads = std::min(max_threads, max_threads_env);
  if (n_max_threads == numbers::invalid_unsigned_int)
    n_max_threads = tbb::task_scheduler_init::default_num_threads();
  else
    {
      static tbb::task_scheduler_init dummy (n_max_threads);
    }
}


unsigned int MultithreadInfo::n_threads() const
{
  if (n_max_threads == numbers::invalid_unsigned_int)
    return tbb::task_scheduler_init::default_num_threads();
  else
    return n_max_threads;
}


#else                            // not in MT mode

unsigned int MultithreadInfo::get_n_cpus()
{
  return 1;
}

unsigned int MultithreadInfo::n_threads() const
{
  return 1;
}

void MultithreadInfo::set_thread_limit(const unsigned int)
{
}

#endif


bool MultithreadInfo::is_running_single_threaded()
{
  return n_threads() == 1;
}


MultithreadInfo::MultithreadInfo ()
  :
  n_cpus (get_n_cpus()),
  n_default_threads (n_cpus),
  n_max_threads (numbers::invalid_unsigned_int)
{}



std::size_t
MultithreadInfo::memory_consumption ()
{
  // only simple data elements, so
  // use sizeof operator
  return sizeof (MultithreadInfo);
}



// definition of the variable which is declared `extern' in the .h file
MultithreadInfo multithread_info;

DEAL_II_NAMESPACE_CLOSE
