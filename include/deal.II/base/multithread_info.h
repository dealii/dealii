// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2015 by the deal.II authors
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

#ifndef dealii__multithread_info_h
#define dealii__multithread_info_h
//---------------------------------------------------------------------------


#include <deal.II/base/config.h>
#include <deal.II/base/types.h>
#include <deal.II/base/exceptions.h>

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
 * explicitly created threads and may want to use a number of threads that is
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
   * The number of CPUs in the system. At the moment detection of CPUs is only
   * implemented on Linux, FreeBSD, and Mac computers.  It is one if detection
   * failed or is not implemented on your system.
   *
   * If it is one, although you are on a multi-processor machine, please refer
   * to the documentation in <tt>multithread_info.cc</tt> near to the
   * <tt>error</tt> directive.
   */
  static unsigned int n_cores ();

  /**
   * @deprecated Use n_cores() instead.
   */
  static const unsigned int n_cpus DEAL_II_DEPRECATED;

  /**
   * Returns the number of threads to use. This is initially set to the number
   * of cores the system has (see n_cores()) but can be further restricted by
   * set_thread_limit().
   */
  static unsigned int n_threads ();

  /**
   * Return an estimate for the memory consumption, in bytes, of this object.
   * This is not exact (but will usually be close) because calculating the
   * memory usage of trees (e.g., <tt>std::map</tt>) is difficult.
   */
  static std::size_t memory_consumption ();

  /**
   * Set the maximum number of threads to be used to the minimum of the
   * environment variable DEAL_II_NUM_THREADS and the given argument (or its
   * default value). This affects the initialization of the TBB. If neither is
   * given, the default from TBB is used (based on the number of cores in the
   * system).
   *
   * Due to limitations in the way TBB can be controlled, only the first call
   * to this method will have any effect. In practice, this means that you
   * need to call this function before you get to any point in your program
   * where multiple threads may be created. In other words, the correct place
   * for a call to this function is at the top of your <code>main()</code>
   * function.
   *
   * This routine is called automatically by MPI_InitFinalize. Use the
   * appropriate argument of the constructor of MPI_InitFinalize if you have
   * an MPI based code.
   */
  static void set_thread_limit (const unsigned int max_threads = numbers::invalid_unsigned_int);


  /**
   * Returns if the TBB is running using a single thread either because of
   * thread affinity or because it is set via a call to set_thread_limit. This
   * is used in the PETScWrappers to avoid using the interface that is not
   * thread-safe.
   */
  static bool is_running_single_threaded ();

  /**
   * Exception
   */
  DeclException0(ExcProcNotPresent);

  /**
   * @deprecated All members are static, so there is no need to construct an
   * instance.
   */
  MultithreadInfo () DEAL_II_DEPRECATED;

private:

  /**
   * Private function to determine the number of CPUs. Implementation for
   * Linux, OSF, SGI, and Sun machines; if no detection of the number of CPUs
   * is supported, or if detection fails, this function returns one.
   */
  static unsigned int get_n_cpus ();

  /**
   * Variable representing the maximum number of threads.
   */
  static unsigned int n_max_threads;
};



/**
 * Global variable of type <tt>MultithreadInfo</tt>.
 *
 * @deprecated Use the static member functions instead.
 *
 * @ingroup threads
 */
extern MultithreadInfo multithread_info DEAL_II_DEPRECATED;




//---------------------------------------------------------------------------
DEAL_II_NAMESPACE_CLOSE
// end of #ifndef dealii__multithread_info_h
#endif
//---------------------------------------------------------------------------
