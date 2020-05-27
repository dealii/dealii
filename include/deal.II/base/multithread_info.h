// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2018 by the deal.II authors
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

#ifndef dealii_multithread_info_h
#  define dealii_multithread_info_h
//---------------------------------------------------------------------------


#  include <deal.II/base/config.h>

#  include <deal.II/base/exceptions.h>
#  include <deal.II/base/types.h>

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
 * @author Thomas Richter, Wolfgang Bangerth, 2000
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

private:
  /**
   * Variable representing the maximum number of threads.
   */
  static unsigned int n_max_threads;
};



//---------------------------------------------------------------------------
DEAL_II_NAMESPACE_CLOSE
// end of #ifndef dealii_multithread_info_h
#endif
//---------------------------------------------------------------------------
