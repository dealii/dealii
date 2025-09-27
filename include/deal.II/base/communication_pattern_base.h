// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_base_communication_pattern_base_h
#define dealii_base_communication_pattern_base_h

#include <deal.II/base/config.h>

#include <deal.II/base/mpi_stub.h>

DEAL_II_NAMESPACE_OPEN

// Forward declaration
#ifndef DOXYGEN
class IndexSet;
#endif

namespace Utilities
{
  namespace MPI
  {
    /**
     * CommunicationPattern is an abstract class that is used to define a
     * communication plan that can be called repeatedly to efficiently obtain
     * off-processor elements. The idea is to decouple the communication pattern
     * from the data that needs to be communicated. The goal is to reuse the
     * same communication pattern for different containers. This is similar to
     * the way SparseMatrix and SparsityPattern works.
     *
     * Conceptually, this class operates under the assumption that data
     * are stored in one linearly indexed array of which every MPI process
     * stores some elements (or possibly none). In practice it does of course
     * not matter whether the elements are stored in contiguous arrays; the
     * point is simply that a single integer index can be used to access a
     * specific element. The elements of this logical array are (or at least
     * may be) stored on different processes in a parallel MPI universe.
     *
     * In this world view, every process has a set of elements and their
     * indices in the array form the "locally owned indices". Next, every
     * process will as part of executing an algorithm require access to some
     * of the elements stored elsewhere; we call the indices of these elements
     * the "ghost indices", in analogy to how vectors and triangulations
     * partition vector elements and mesh cells into locally
     * owned ones and ghost ones (along, of course, with cells and ghosts stored
     * on other processes that the current process simply does not care about
     * and consequently needs not know anything about).
     *
     * The point of this class (and its implementations in derived classes) is
     * to set up communication infrastructure whereby one process can inquire
     * efficiently about the "ghost elements" stored on other processes, and
     * to send those locally owned elements to those other processes that
     * require knowledge of their value because they list these elements among
     * their respective "ghost indices".
     */
    class CommunicationPatternBase
    {
    public:
      /**
       * Destructor.
       */
      virtual ~CommunicationPatternBase() = default;

      /**
       * Reinitialize the communication pattern.
       *
       * @param[in] locally_owned_indices The set of indices of elements
       *   in the array mentioned in the class documentation that are
       *   stored on the current process.
       * @param[in] ghost_indices The set of indices of elements in the
       *   array mentioned in the class documentation that the current
       *   process will need to be able to import.
       * @param[in] communicator The MPI communicator used to describe the
       *   entire set of processes that participate in the storage and
       *   access to elements of the array.
       */
      virtual void
      reinit(const IndexSet &locally_owned_indices,
             const IndexSet &ghost_indices,
             const MPI_Comm  communicator) = 0;

      /**
       * Return the underlying MPI communicator.
       */
      virtual MPI_Comm
      get_mpi_communicator() const = 0;
    };

  } // namespace MPI

} // namespace Utilities

DEAL_II_NAMESPACE_CLOSE

#endif
