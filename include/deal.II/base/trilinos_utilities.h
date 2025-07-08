// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_trilinos_utilities_h
#define dealii_trilinos_utilities_h

#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>

#ifdef DEAL_II_WITH_TRILINOS

DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
#  include <Epetra_Comm.h>
#  include <Epetra_Map.h>
#  include <Teuchos_Comm.hpp>
#  include <Teuchos_RCP.hpp>

#  ifdef DEAL_II_WITH_MPI
#    include <Epetra_MpiComm.h>
#  else
#    include <Epetra_SerialComm.h>
#  endif

#  ifdef DEAL_II_TRILINOS_WITH_TPETRA
#    include <Teuchos_RCPDecl.hpp>
#  endif // DEAL_II_TRILINOS_WITH_TPETRA
DEAL_II_ENABLE_EXTRA_DIAGNOSTICS
#endif


DEAL_II_NAMESPACE_OPEN

/**
 * A namespace for utility functions that are not particularly specific to
 * finite element computing or numerical programs, but nevertheless are needed
 * in various contexts when writing applications.
 *
 * @ingroup utilities
 */
namespace Utilities
{
#ifdef DEAL_II_WITH_TRILINOS
  /**
   * This namespace provides some of the basic structures used in the
   * initialization of the Trilinos objects (e.g., matrices, vectors, and
   * preconditioners).
   */
  namespace Trilinos
  {
    /**
     * Return a Trilinos Epetra_Comm object needed for creation of
     * Epetra_Maps.
     *
     * If deal.II has been configured to use a compiler that does not support
     * MPI then the resulting communicator will be a serial one. Otherwise,
     * the communicator will correspond to MPI_COMM_WORLD, i.e. a communicator
     * that encompasses all processes within this MPI universe.
     */
    const Epetra_Comm &
    comm_world();

    /**
     * Return a Trilinos Epetra_Comm object needed for creation of
     * Epetra_Maps.
     *
     * If deal.II has been configured to use a compiler that does not support
     * MPI then the resulting communicator will be a serial one. Otherwise,
     * the communicator will correspond to MPI_COMM_SELF, i.e. a communicator
     * that comprises only this one processor.
     */
    const Epetra_Comm &
    comm_self();

    /**
     * Return a Teuchos::Comm object needed for creation of Tpetra::Maps.
     *
     * If deal.II has been configured to use a compiler that does not support
     * MPI then the resulting communicator will be a serial one. Otherwise,
     * the communicator will correspond to MPI_COMM_SELF, i.e. a communicator
     * that comprises only this one processor.
     */
    const Teuchos::RCP<const Teuchos::Comm<int>> &
    tpetra_comm_self();

    /**
     * Given a communicator, duplicate it. If the given communicator is
     * serial, that means to just return a copy of itself. On the other hand,
     * if it is %parallel, we duplicate the underlying MPI_Comm object: we
     * create a separate MPI communicator that contains the same processors
     * and in the same order but has a separate identifier distinct from the
     * given communicator. The function returns a pointer to a new object of a
     * class derived from Epetra_Comm. The caller of this function needs to
     * assume ownership of this function. The returned object should be
     * destroyed using the destroy_communicator() function.
     *
     * This facility is used to separate streams of communication. For
     * example, a program could simply use MPI_Comm_World for everything. But
     * it is easy to come up with scenarios where sometimes not all processors
     * participate in a communication that is intended to be global -- for
     * example if we assemble a matrix on a coarse mesh with fewer cells than
     * there are processors, some processors may not sync their matrices with
     * the rest because they haven't written into it because they own no
     * cells. That's clearly a bug. However, if these processors just continue
     * their work, and the next %parallel operation happens to be a sync on a
     * different matrix, then the sync could succeed -- by accident, since
     * different processors are talking about different matrices.
     *
     * This kind of situation can be avoided if we use different communicators
     * for different matrices which reduces the likelihood that communications
     * meant to be separate aren't recognized as such just because they happen
     * on the same communicator. In addition, it is conceivable that some MPI
     * operations can be parallelized using multiple threads because their
     * communicators identifies the communication in question, not their
     * relative timing as is the case in a sequential program that just uses a
     * single communicator.
     */
    Epetra_Comm *
    duplicate_communicator(const Epetra_Comm &communicator);

    /**
     * Given an Epetra communicator that was created by the
     * duplicate_communicator() function, destroy the underlying MPI
     * communicator object and reset the Epetra_Comm object to a the result of
     * comm_self().
     *
     * It is necessary to call this function at the time when the result of
     * duplicate_communicator() is no longer needed. The reason is that in
     * that function, we first create a new MPI_Comm object and then create an
     * Epetra_Comm around it. While we can take care of destroying the latter,
     * it doesn't destroy the communicator since it can only assume that it
     * may also be still used by other objects in the program. Consequently,
     * we have to take care of destroying it ourselves, explicitly.
     *
     * This function does exactly that. Because this has to happen while the
     * Epetra_Comm object is still around, it first resets the latter and then
     * destroys the communicator object.
     *
     * @note If you call this function on an Epetra_Comm object that is not
     * created by duplicate_communicator(), you are likely doing something
     * quite wrong. Don't do this.
     */
    void
    destroy_communicator(Epetra_Comm &communicator);

    /**
     * Return the number of MPI processes there exist in the given
     * @ref GlossMPICommunicator "communicator"
     * object. If this is a sequential job (i.e., the program
     * is not using MPI at all, or is using MPI but has been started with
     * only one MPI process), then the communicator necessarily involves
     * only one process and the function returns 1.
     */
    unsigned int
    get_n_mpi_processes(const Epetra_Comm &mpi_communicator);

    /**
     * Return the number of the present MPI process in the space of processes
     * described by the given communicator. This will be a unique value for
     * each process between zero and (less than) the number of all processes
     * (given by get_n_mpi_processes()).
     */
    unsigned int
    get_this_mpi_process(const Epetra_Comm &mpi_communicator);

    /**
     * Given a Trilinos Epetra map, create a new map that has the same
     * subdivision of elements to processors but uses the given communicator
     * object instead of the one stored in the first argument. In essence,
     * this means that we create a map that communicates among the same
     * processors in the same way, but using a separate channel.
     *
     * This function is typically used with a communicator that has been
     * obtained by the duplicate_communicator() function.
     */
    Epetra_Map
    duplicate_map(const Epetra_BlockMap &map, const Epetra_Comm &comm);
  } // namespace Trilinos
#endif


#ifdef DEAL_II_TRILINOS_WITH_TPETRA
  namespace Trilinos
  {
    /**
     * Return the underlying MPI_Comm communicator from the
     * <a
     * href="https://docs.trilinos.org/dev/packages/teuchos/doc/html/classTeuchos_1_1Comm.html">Teuchos::Comm</a>
     * communicator.
     */
    MPI_Comm
    teuchos_comm_to_mpi_comm(
      const Teuchos::RCP<const Teuchos::Comm<int>> &teuchos_comm);

    namespace internal
    {
      /**
       * Creates and returns a
       * <a
       * href="https://docs.trilinos.org/dev/packages/teuchos/doc/html/classTeuchos_1_1RCP.html">Teuchos::RCP</a>
       * object for type T.
       *
       * @note In Trilinos 14.0.0, the function
       * <a
       * href="https://docs.trilinos.org/dev/packages/teuchos/doc/html/namespaceTeuchos.html#a280c0ab8c9ee8d0481114d4edf5a3393">Teuchos::make_rcp()</a>
       * was introduced, which should be preferred to this function.
       */
#  if defined(DOXYGEN) || !DEAL_II_TRILINOS_VERSION_GTE(14, 0, 0)
      template <class T, class... Args>
      Teuchos::RCP<T>
      make_rcp(Args &&...args);
#  else
      using Teuchos::make_rcp;
#  endif // defined DOXYGEN || !DEAL_II_TRILINOS_VERSION_GTE(14, 0, 0)
    }    // namespace internal



    /* ------------------------- Inline functions ---------------------- */
    namespace internal
    {
#  if !DEAL_II_TRILINOS_VERSION_GTE(14, 0, 0)
      template <class T, class... Args>
      inline Teuchos::RCP<T>
      make_rcp(Args &&...args)
      {
        return Teuchos::RCP<T>(new T(std::forward<Args>(args)...));
      }
#  endif // !DEAL_II_TRILINOS_VERSION_GTE(14, 0, 0)
    }    // namespace internal

  }    // namespace Trilinos
#endif // DEAL_II_TRILINOS_WITH_TPETRA

} // namespace Utilities

DEAL_II_NAMESPACE_CLOSE

#endif
