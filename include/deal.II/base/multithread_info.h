// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2000 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_multithread_info_h
#  define dealii_multithread_info_h
//---------------------------------------------------------------------------


#  include <deal.II/base/config.h>

#  include <deal.II/base/exceptions.h>
#  include <deal.II/base/types.h>

#  include <memory>

#  ifdef DEAL_II_WITH_TASKFLOW
// forward declaration from <taskflow/taskflow.hpp>
namespace tf
{
  class Executor;
}
#  endif


DEAL_II_NAMESPACE_OPEN

/**
 * This class provides information about the system which may be of use in
 * multithreaded programs.  At the moment this is just the number of CPUs. If
 * deal.II is compiled with multithreading support, some functions will use
 * multiple threads for their action. Currently the library supports both
 * thread-based and task-based parallelism.
 * @ref threads
 * describes the different uses of each. The default number of threads used
 * for task-based parallel methods is selected automatically by the Threading
 * Building Blocks library. See
 * @ref threads
 * for more information on this.  Thread-based parallel methods need to
 * explicitly create threads and may want to use a number of threads that is
 * related to the number of CPUs in your system. The recommended number of
 * threads can be queried using MultithreadInfo::n_threads(), while the number
 * of cores in the system is returned by MultithreadInfo::n_cores().
 *
 * @ingroup threads
 */
class MultithreadInfo
{
public:
  /**
   * Constructor. This constructor is deleted because no instance of
   * this class needs to be constructed (all members are static).
   */
  MultithreadInfo() = delete;

  /**
   * The number of CPUs in the system.
   *
   * This internally calls
   * [<code>std::thread::hardware_concurrency</code>](https://en.cppreference.com/w/cpp/thread/thread/hardware_concurrency)
   * but sets the result to 1 if the call returns an error.
   */
  static unsigned int
  n_cores();

  /**
   * Return the number of threads to use. This is initially set to the number
   * of cores the system has (see n_cores()) but can be further restricted by
   * set_thread_limit() and the environment variable DEAL_II_NUM_THREADS.
   */
  static unsigned int
  n_threads();

  /**
   * Return an estimate for the memory consumption, in bytes, of this object.
   * This is not exact (but will usually be close) because calculating the
   * memory usage of trees (e.g., <tt>std::map</tt>) is difficult.
   */
  static std::size_t
  memory_consumption();

  /**
   * Set the maximum number of threads to be used to the minimum of the
   * environment variable DEAL_II_NUM_THREADS and the given argument (or its
   * default value). This affects the initialization of the TBB. If neither is
   * given, the default from TBB is used (based on the number of cores in the
   * system).
   *
   * This routine is executed automatically with the default argument before
   * your code in main() is running (using a static constructor). It is also
   * executed by Utilities::MPI::MPI_InitFinalize. Use the appropriate
   * argument of the constructor of Utilities::MPI::MPI_InitFinalize if you
   * have an MPI based code.
   */
  static void
  set_thread_limit(
    const unsigned int max_threads = numbers::invalid_unsigned_int);

  /**
   * Return if the TBB is running using a single thread either because of
   * thread affinity or because it is set via a call to set_thread_limit. This
   * is used in the PETScWrappers to avoid using the interface that is not
   * thread-safe.
   */
  static bool
  is_running_single_threaded();

  /**
   * Make sure the multithreading API is initialized. This normally does not
   * need to be called in usercode.
   */
  static void
  initialize_multithreading();


#  ifdef DEAL_II_WITH_TASKFLOW
  /**
   * Return a reference to the global Executor from taskflow.
   *
   * The Executor is set to use n_threads() worker threads that you can
   * control using set_thread_limit() and the DEAL_II_NUM_THREADS environment
   * variable.
   */
  static tf::Executor &
  get_taskflow_executor();
#  endif

private:
  /**
   * Variable representing the maximum number of threads.
   */
  static unsigned int n_max_threads;

#  ifdef DEAL_II_WITH_TASKFLOW
  /**
   * Store a taskflow Executor that is constructed with N workers (from
   * set_thread_limit).
   */
  static std::unique_ptr<tf::Executor> executor;
#  endif
};



//---------------------------------------------------------------------------
DEAL_II_NAMESPACE_CLOSE
// end of #ifndef dealii_multithread_info_h
#endif
//---------------------------------------------------------------------------
