// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_init_finalize_h
#define dealii_init_finalize_h

#include <deal.II/base/config.h>

#include <deal.II/base/mpi_stub.h>
#include <deal.II/base/numbers.h>
#include <deal.II/base/types.h>

#include <boost/signals2.hpp>

#include <set>

DEAL_II_NAMESPACE_OPEN

/**
 * The enum type given to the constructor of InitFinalize telling which
 * external libraries should be initialized and finalized by deal.II.
 */
enum class InitializeLibrary
{
  /**
   * No external library will be initialized/finalized.
   */
  None = 0,
  /**
   * Initialize/finalize MPI.
   */
  MPI = 1,
  /**
   * Initialize/finalize Kokkos.
   */
  Kokkos = 2,
  /**
   * Initialize/finalize both PETSc and SLEPc.
   */
  SLEPc = 4,
  /**
   * Initialize/finalize PETSc.
   */
  PETSc = 8,
  /**
   * Initialize/finalize Zoltan.
   */
  Zoltan = 16,
  /**
   * Initialize/finalize P4EST and SC.
   */
  P4EST = 32
};



/**
 * Global operator which returns an object in which all bits are set which are
 * either set in the first or the second argument.
 */
inline InitializeLibrary
operator|(const InitializeLibrary f1, const InitializeLibrary f2)
{
  return static_cast<InitializeLibrary>(static_cast<unsigned int>(f1) |
                                        static_cast<unsigned int>(f2));
}



/**
 * Global operator which returns an object in which all bits are set that are
 * set both in the first and second argument.
 */
inline InitializeLibrary
operator&(const InitializeLibrary f1, const InitializeLibrary f2)
{
  return static_cast<InitializeLibrary>(static_cast<unsigned int>(f1) &
                                        static_cast<unsigned int>(f2));
}



/**
 * A class that is used to initialize and finalize deal.II and the external
 * libraries requested by the user.
 */
class InitFinalize
{
public:
  /**
   * Initialize deal.II and the requested external libraries.
   *
   * @note This function calls MultithreadInfo::set_thread_limit() with
   * either @p max_num_threads or, following the discussion above, a
   * number of threads equal to the number of cores allocated to this MPI
   * process. However, MultithreadInfo::set_thread_limit() in turn also
   * evaluates the environment variable DEAL_II_NUM_THREADS. Finally, the
   * worker threads can only be created on cores to which the current MPI
   * process has access to; some MPI implementations limit the number of
   * cores each process may access to one or a subset of cores in order to
   * ensure better cache behavior. Consequently, the number of threads
   * that will really be created will be the minimum of the argument
   * passed here, the environment variable (if set), and the number of
   * cores accessible to the thread.
   *
   * @note MultithreadInfo::set_thread_limit() can only work if it is
   * called before any threads are created. The safest place for a call to
   * it is therefore at the beginning of <code>main()</code>.
   * Consequently, this extends to the current class: the best place to
   * create an object of this type is also at or close to the top of
   * <code>main()</code>.
   */
  InitFinalize(
    int                     &argc,
    char                  **&argv,
    const InitializeLibrary &libraries,
    const unsigned int       max_num_threads = numbers::invalid_unsigned_int);

  /**
   * Destructor. Finalize deal.II and the external libraries initialize by
   * deal.II.
   */
  ~InitFinalize();

  /**
   * Finalize deal.II and the external libraries initialize by
   * deal.II.
   */
  void
  finalize();

  /**
   * Register a reference to an MPI_Request
   * on which we need to call `MPI_Wait` before calling `MPI_Finalize`.
   *
   * The object @p request needs to exist when MPI_Finalize is called, which means the
   * request is typically statically allocated. Otherwise, you need to call
   * unregister_request() before the request goes out of scope. Note that it
   * is acceptable for a request to be already waited on (and consequently
   * reset to MPI_REQUEST_NULL).
   *
   * It is acceptable to call this function more than once with the same
   * instance (as it is done in the example below).
   *
   * Typically, this function is used by CollectiveMutex and not directly,
   * but it can also be used directly like this:
   * @code
   * void my_fancy_communication()
   * {
   *   static MPI_Request request = MPI_REQUEST_NULL;
   *   MPI_InitFinalize::register_request(request);
   *   MPI_Wait(&request, MPI_STATUS_IGNORE);
   *   // [some algorithm that is not safe to be executed twice in a row.]
   *   MPI_IBarrier(comm, &request);
   * }
   * @endcode
   */
  static void
  register_request(MPI_Request &request);

  /**
   * Unregister a request previously added using register_request().
   */
  static void
  unregister_request(MPI_Request &request);

  /**
   * A structure that has boost::signal objects to register a call back
   * to run after MPI init or finalize.
   *
   * For documentation on signals, see
   * http://www.boost.org/doc/libs/release/libs/signals2 .
   */
  struct Signals
  {
    /**
     * A signal that is triggered immediately after we have
     * initialized the MPI context with <code>MPI_Init()</code>.
     */
    boost::signals2::signal<void()> at_mpi_init;

    /**
     * A signal that is triggered just before we close the MPI context
     * with <code>MPI_Finalize()</code>. It can be used to deallocate
     * statically allocated MPI resources that need to be deallocated
     * before <code>MPI_Finalize()</code> is called.
     */
    boost::signals2::signal<void()> at_mpi_finalize;
  };

  static Signals signals;

private:
  InitializeLibrary libraries;

  /**
   * Requests to MPI_Wait before finalizing
   */
  static std::set<MPI_Request *> requests;

  bool is_finalized = false;

#ifdef DEAL_II_WITH_PETSC
  bool finalize_petscslepc;
#endif
};

DEAL_II_NAMESPACE_CLOSE

#endif
