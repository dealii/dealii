// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2012 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_partitioner_h
#define dealii_partitioner_h

#include <deal.II/base/config.h>

#include <deal.II/base/array_view.h>
#include <deal.II/base/communication_pattern_base.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/memory_space.h>
#include <deal.II/base/mpi_stub.h>
#include <deal.II/base/types.h>

#include <deal.II/lac/vector_operation.h>

#include <Kokkos_Core.hpp>

#include <limits>
#include <memory>
#include <vector>


DEAL_II_NAMESPACE_OPEN

namespace Utilities
{
  namespace MPI
  {
    /**
     * This class defines a model for the partitioning of a vector (or, in
     * fact, any linear data structure) among processors using MPI.
     *
     * The partitioner stores the global vector size and the locally owned
     * range as a half-open interval [@p lower, @p upper) on each process. In
     * other words, it assumes that every process stores a contiguous subset
     * of the array mentioned in the documentation of the base class
     * Utilities::MPI::CommunicationPatternBase. (If you want to store
     * non-contiguous parts of these arrays on each process, take a look
     * at Utilities::MPI::NoncontiguousPartitioner.)
     *
     * Furthermore, it includes a structure for the point-to-point communication
     * patterns. It allows the inclusion of ghost indices (i.e. indices that a
     * current processor needs to have access to, but are owned by another
     * process) through an IndexSet. In addition, it also stores the other
     * processors' ghost indices belonging to the current processor (see
     * import_targets()), which are the indices where other processors might
     * require information from. In a sense, these import indices form the dual
     * of the ghost indices. This information is gathered once when constructing
     * the partitioner, which obviates subsequent global communication steps
     * when exchanging data.
     *
     * The figure below gives an example of index space $[0,74)$ being split
     * into four parts that are each owned by one MPI process:
     * @image html partitioner.png
     * The first row (above the thick black line) shows which process owns which
     * elements. Below it, the next four lines indicate which elements of
     * the overall array each processor wants to know about -- this is
     * generally a superset of the locally owned elements, with the difference
     * being what are called "ghost elements".
     *
     * To understand the remaining pieces of the figure (and this class),
     * remember that in MPI, you can't just ask another process for data.
     * (That's not quite true: There are mechanisms in newer MPI standards
     * for this, but as a general rule it's true.) Rather, if you need
     * information, you need to send another process a message, the other
     * process needs to expect the message and respond as appropriate with
     * a message of its own. In practice, it is therefore easier and faster
     * if each process will *already* know what it will be asked and, at
     * the appropriate time, just send that data. The remaining lines of
     * information set up this kind of scheme.
     *
     * To this end, note that process 0 will want to know about five
     * elements it does not own itself: those with indices 20, 21 (owned
     * by process 1); and 40, 41, 43 (owned by process 2). Similar information
     * can be obtained from the following lines. To satisfy this need for
     * knowledge, it would therefore be quite useful if process 1 stored
     * that, at the appropriate time, it will have to send elements 20, 21
     * to process 0. Similarly, if you go through lines 2-4 below the thick
     * black line, you will see that process 0 should know that it will need
     * to send elements 1, 2, 13, 18, 19 to process 1; 18, 19 to process 2;
     * and elements 1, 2, 13 to process 3. These are called "import indices"
     * because other processes will want to import them. Instead of storing
     * these indices as a set, it is often useful to use half-open index
     * sets instead, and so the import indices listed above form the following
     * collection of sets: `[1,3)`, `[13,14)`, `[18,20)`, `[18,20)`,
     * `[1,3)`, `[13,14)`. This is how the import indices are shown
     * in the figure above. We now only have to know which of these
     * half-open sets are to be sent to whom. This is done in the line
     * above it where we list the "import targets" (i.e., the target
     * processes for an import operations): Process 1 will receive 5
     * elements (which are comprised of the first three half-open
     * target index sets), process 2 will receive 2 indices
     * (the fourth half-open interval), and process 3 will receive
     * 3 indices (the remaining two half-open intervals). This information
     * is encoded as the pairs `{1,5}`, `{2,2}`, `{3,3}` as the import
     * targets. Similar considerations can be made for what processes
     * 1, 2, and 3 will have to send out.
     *
     * Finally, when receiving information, it is useful to know how
     * many indices each process will receive from whom since then one
     * can already pre-allocate buffers of the right size. This is listed
     * in the last line under "ghost targets": Process 0 will receive
     * two elements from process 1 (namely those with indices 20, 21),
     * and three from process 2 (namely those with indices 40, 41, 43).
     * This is encoded as pairs `{1,2}` and `{2,3}`. Again, similar
     * considerations can be made for what processes 1, 2, and 3
     * should expect, and what is then shown in their respective columns.
     *
     * The main purpose of this class is to set up these data structures
     * knowing only which process owns which elements, and for which
     * additional ghost elements each process needs knowledge.
     *
     *
     * <h4>Local and global numbering</h4>
     *
     * The partitioner includes a mechanism for converting global to local and
     * local to global indices. Internally, this class stores vector elements
     * using the convention as follows: The local range is associated with
     * local indices [0, locally_owned_size()), and ghost indices are stored
     * consecutively in [locally_owned_size(), locally_owned_size() +
     * n_ghost_indices()). The ghost indices are sorted according to their
     * global index.
     *
     *
     * <h4>Parallel data exchange</h4>
     *
     * This class also handles the ghost data exchange for partitioned
     * arrays of objects -- i.e., if a separate class stores the
     * locally owned elements on each process, then this class
     * facilitates the importation of those elements that are locally
     * needed as ghosts but stored elsewhere. An example of where this class
     * is used is the LinearAlgebra::distributed::Vector class.
     *
     * The data exchange happens through four functions:
     * <ul>
     * <li> export_to_ghosted_array_start() is used for initiating an export
     * operation that sends data from the locally owned data field, passed as
     * an array, to the ghost data arrays of other processes (according to the
     * ghost indices stored in the present class).
     * This call starts non-blocking MPI send communication
     * routines, but does not wait for the routines to finish. Thus, the user
     * may not write into the respective positions of the underlying arrays as
     * the data might still be needed by MPI.
     * <li> export_to_ghosted_array_finish() finalizes the MPI data exchange
     * started in export_to_ghosted_array_start() and signals that the data in
     * the arrays may be used for further processing or modified as
     * appropriate.
     * <li> import_from_ghosted_array_start() is used for initiating an import
     * operation that sends data from a ghost data field, passed as an array,
     * to the locally owned array according to the ghost indices stored in the
     * present class. A VectorOperation::values flag can be passed to decide
     * on how to combine the data in the ghost field with the data at the
     * owner, since both relate to the same data entry. In assembly, this is
     * usually an add-to operation where contributions from all processes to
     * a locally owned element need to be added up. This call starts
     * non-blocking MPI communication routines, but does not wait for the
     * routines to finish. Thus, the user may not write into the respective
     * positions of the underlying arrays as the data might still be needed by
     * MPI.
     * <li> import_from_ghosted_array_finish() finalizes the MPI data exchange
     * started in import_from_ghosted_array_start() and signals that the data
     * in the arrays may be used for further processing or modified as
     * appropriate.
     * </ul>
     *
     * The MPI communication routines are point-to-point communication patterns.
     *
     *
     * <h4>Sending only selected ghost data</h4>
     *
     * This partitioner class operates on a fixed set of ghost indices and
     * must always be compatible with the ghost indices inside the array whose
     * partitioning it represents. In some cases,
     * one only wants to send around some of the ghost indices present in
     * a vector, but without creating a copy of the vector with a
     * suitable index set - think e.g. of local time stepping where different
     * regions of a vector might be exchanged at different stages of a time
     * step slice. This class supports that case by the following model: A
     * vector is first created with the full ghosted index set. Then, a second
     * Partitioner instance is created that sets ghost indices with a tighter
     * index set as ghosts, but specifying the larger index set as the second
     * argument to the set_ghost_indices() call. When data is exchanged, the
     * export_to_ghosted_array_start() and import_from_ghosted_array_start()
     * detect this case and only send the selected indices, taken from the
     * full array of ghost entries.
     */
    class Partitioner : public Utilities::MPI::CommunicationPatternBase
    {
    public:
      /**
       * Default constructor.
       */
      Partitioner();

      /**
       * Constructor with size argument. Creates an MPI_COMM_SELF structure
       * where there is no real parallel layout.
       */
      Partitioner(const unsigned int size);

      /**
       * Constructor that takes the number of locally-owned degrees of freedom
       * @p local_size and the number of ghost degrees of freedom @p ghost_size.
       *
       * The local index range is translated to global indices in an ascending
       * and one-to-one fashion, i.e., the indices of process $p$ sit exactly
       * between the indices of the processes $p-1$ and $p+1$, respectively.
       *
       * @note Setting the @p ghost_size variable to an appropriate value
       *   provides memory space for the ghost data in a vector's memory
       *   allocation as and allows access to it via local_element(). However,
       *   the associated global indices must be handled externally in this
       *   case.
       */
      Partitioner(const types::global_dof_index local_size,
                  const types::global_dof_index ghost_size,
                  const MPI_Comm                communicator);

      /**
       * Constructor with index set arguments. This constructor creates a
       * distributed layout based on a given communicators, an IndexSet
       * describing the locally owned range and another one for describing
       * ghost indices that are owned by other processors, but that we need to
       * have read or write access to.
       */
      Partitioner(const IndexSet &locally_owned_indices,
                  const IndexSet &ghost_indices_in,
                  const MPI_Comm  communicator_in);

      /**
       * Constructor with one index set argument. This constructor creates a
       * distributed layout based on a given communicator, and an IndexSet
       * describing the locally owned range. It allows to set the ghost
       * indices at a later time. Apart from this, it is similar to the other
       * constructor with two index sets.
       */
      Partitioner(const IndexSet &locally_owned_indices,
                  const MPI_Comm  communicator_in);

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
             const MPI_Comm  communicator) override;

      /**
       * Set the locally owned indices. Used in the constructor.
       */
      void
      set_owned_indices(const IndexSet &locally_owned_indices);

      /**
       * Set the ghost indices after the constructor has been
       * called.
       *
       * The optional parameter @p larger_ghost_index_set allows defining an
       * indirect addressing into a larger set of ghost indices. This setup is
       * useful if a distributed vector is based on that larger ghost index
       * set but only a tighter subset should be communicated according to
       * @p ghost_indices.
       */
      void
      set_ghost_indices(const IndexSet &ghost_indices,
                        const IndexSet &larger_ghost_index_set = IndexSet());

      /**
       * Return the global size.
       */
      types::global_dof_index
      size() const;

      /**
       * Return the number of locally owned indices,
       * i.e., local_range().second minus local_range().first.
       * The returned numbers need to add up to the total number of indices when
       * summed over all processes
       */
      unsigned int
      locally_owned_size() const;

      /**
       * Return an IndexSet representation of the local range. This class
       * only supports contiguous local ranges, so the IndexSet actually only
       * consists of one single range of data, and is equivalent to the result
       * of local_range().
       */
      const IndexSet &
      locally_owned_range() const;

      /**
       * Return the local range. The returned pair consists of the index of
       * the first element and the index of the element one past the last
       * locally owned one.
       */
      std::pair<types::global_dof_index, types::global_dof_index>
      local_range() const;

      /**
       * Return true if the given global index is in the local range of this
       * processor.
       */
      bool
      in_local_range(const types::global_dof_index global_index) const;

      /**
       * Return the local index corresponding to the given global index. If
       * the given global index is neither locally owned nor a ghost, an
       * exception is thrown.
       *
       * Note that the returned local index for locally owned indices will be
       * between 0 and locally_owned_size() - 1, and the local index for
       * ghosts is between locally_owned_size() and locally_owned_size() +
       * n_ghost_indices() - 1.
       */
      unsigned int
      global_to_local(const types::global_dof_index global_index) const;

      /**
       * Return the global index corresponding to the given local index.
       *
       * Note that the local index for locally owned indices must be between 0
       * and locally_owned_size() - 1, and the local index for ghosts must be
       * between locally_owned_size() and locally_owned_size() +
       * n_ghost_indices() - 1.
       */
      types::global_dof_index
      local_to_global(const unsigned int local_index) const;

      /**
       * Return whether the given global index is a ghost index on the
       * present processor. Returns false for indices that are owned locally
       * and for indices not present at all.
       */
      bool
      is_ghost_entry(const types::global_dof_index global_index) const;

      /**
       * Return an IndexSet representation of all ghost indices.
       */
      const IndexSet &
      ghost_indices() const;

      /**
       * Return the number of ghost indices. Same as
       * ghost_indices().n_elements(), but cached for simpler access.
       */
      unsigned int
      n_ghost_indices() const;

      /**
       * In case the partitioner was built to define ghost indices as a subset
       * of indices in a larger set of ghosts, this function returns the
       * numbering in terms of ranges within that set. Similar structure as in
       * an IndexSet, but tailored to be iterated over.
       *
       * In case the partitioner did not take a second set of ghost indices
       * into account, this subset is simply defined as the half-open interval
       * <code>[0, n_ghost_indices())</code>.
       */
      const std::vector<std::pair<unsigned int, unsigned int>> &
      ghost_indices_within_larger_ghost_set() const;

      /**
       * Return a list of processors (first entry) and the number of ghost
       * degrees of freedom owned by that processor (second entry). The sum of
       * the latter over all processors equals n_ghost_indices().
       */
      const std::vector<std::pair<unsigned int, unsigned int>> &
      ghost_targets() const;

      /**
       * Return a vector of ranges of local indices that we are importing during
       * compress(), i.e., others' ghosts that belong to the local range.
       * Similar structure as in an IndexSet, but tailored to be iterated over,
       * and some indices may be duplicated. The returned pairs consists of the
       * index of the first element and the index of the element one past the
       * last one in a range.
       */
      const std::vector<std::pair<unsigned int, unsigned int>> &
      import_indices() const;

      /**
       * Number of import indices, i.e., indices that are ghosts on other
       * processors and we will receive data from.
       */
      unsigned int
      n_import_indices() const;

      /**
       * Return a list of processors (first entry) and the number of degrees
       * of freedom imported from it during compress() operation (second entry)
       * for all the processors that data is obtained from, i.e., locally owned
       * indices that are ghosts on other processors.
       *
       * @note The returned vector only contains those processor id's for which
       * the second entry is non-zero.
       */
      const std::vector<std::pair<unsigned int, unsigned int>> &
      import_targets() const;

      /**
       * Check whether the given partitioner is compatible with the current
       * partitioner. Two partitioners are compatible if
       * they have the same local sizes and the same ghost indices. They do not
       * necessarily need to correspond to the same data that is stored based
       * on these partioner objects. This is a local operation
       * only, i.e., if only some processors decide that the partitioning is
       * not compatible, only these processors will return @p false, whereas
       * the other processors will return @p true.
       */
      bool
      is_compatible(const Partitioner &part) const;

      /**
       * Check whether the given partitioner is compatible with the
       * current partitioner. Two partitioners are compatible if
       * they have the same local size and the same ghost indices. They do not
       * necessarily need to correspond to the same data that is stored based
       * on these partioner objects. As opposed to
       * is_compatible(), this method checks for compatibility among all
       * processors and the method only returns @p true if the partitioner is
       * the same on all processors. In other words, it does a global
       * "and" operation over the results returned by is_compatible() on all
       * involved processes.
       *
       * This method performs global communication, so make sure to use it
       * only in a context where all processors call it the same number of
       * times.
       */
      bool
      is_globally_compatible(const Partitioner &part) const;

      /**
       * Return the MPI ID of the calling processor. Cached to have simple
       * access.
       */
      unsigned int
      this_mpi_process() const;

      /**
       * Return the total number of MPI processor participating in the given
       * partitioner. Cached to have simple access.
       */
      unsigned int
      n_mpi_processes() const;

      /**
       * Return the underlying MPI communicator.
       */
      virtual MPI_Comm
      get_mpi_communicator() const override;

      /**
       * Return whether ghost indices have been explicitly added as a @p
       * ghost_indices argument. Only true if a reinit() call or constructor
       * provided that argument.
       */
      bool
      ghost_indices_initialized() const;

#ifdef DEAL_II_WITH_MPI
      /**
       * Start the exportation of the data in a locally owned array to the
       * range described by the ghost indices of this class.
       *
       * @param communication_channel Sets an offset to the MPI_Isend and
       * MPI_Irecv calls that avoids interference with other ongoing
       * export_to_ghosted_array_start() calls on different entries. Typically
       * handled within the blocks of a block vector. Any value less than 200
       * is a valid value.
       *
       * @param locally_owned_array The array of data from which the data is
       * extracted and sent to the ghost entries on a remote processor.
       *
       * @param temporary_storage A temporary storage array of length
       * n_import_indices() that is used to hold the packed data from the @p
       * locally_owned_array to be sent. Note that this array must not be
       * touched until the respective export_to_ghosted_array_finish() call
       * has been made because the model uses non-blocking communication.
       *
       * @param ghost_array The array that will receive the exported data,
       * i.e., the entries that a remote processor sent to the calling
       * process. Its size must either be n_ghost_indices() or equal the
       * number of ghost indices in the larger index set that was given as
       * second argument to set_ghost_indices(). In case only selected indices
       * are sent, no guarantee is made regarding the entries that do not get
       * set. Some of them might be used to organize the transfer and later
       * reset to zero, so make sure you do not use them in computations.
       *
       * @param requests The list of MPI requests for the ongoing non-blocking
       * communication that will be finalized in the
       * export_to_ghosted_array_finish() call.
       *
       * This functionality is used in
       * LinearAlgebra::distributed::Vector::update_ghost_values().
       */
      template <typename Number, typename MemorySpaceType = MemorySpace::Host>
      void
      export_to_ghosted_array_start(
        const unsigned int                              communication_channel,
        const ArrayView<const Number, MemorySpaceType> &locally_owned_array,
        const ArrayView<Number, MemorySpaceType>       &temporary_storage,
        const ArrayView<Number, MemorySpaceType>       &ghost_array,
        std::vector<MPI_Request>                       &requests) const;

      /**
       * Finish the exportation of the data in a locally owned array to the
       * range described by the ghost indices of this class.
       *
       * @param ghost_array The array that will receive the exported data
       * started in the @p export_to_ghosted_array_start(). This must be the
       * same array as passed to that function, otherwise the behavior is
       * undefined.
       *
       * @param requests The list of MPI requests for the ongoing non-blocking
       * communication that were started in the
       * export_to_ghosted_array_start() call. This must be the same array as
       * passed to that function, otherwise MPI will likely throw an error.
       *
       * This functionality is used in
       * LinearAlgebra::distributed::Vector::update_ghost_values().
       */
      template <typename Number, typename MemorySpaceType = MemorySpace::Host>
      void
      export_to_ghosted_array_finish(
        const ArrayView<Number, MemorySpaceType> &ghost_array,
        std::vector<MPI_Request>                 &requests) const;

      /**
       * Start importing the data on an array indexed by the ghost indices of
       * this class that is later accumulated into a locally owned array with
       * import_from_ghosted_array_finish().
       *
       * @param vector_operation Defines how the data sent to the owner should
       * be combined with the existing entries, e.g., added into.
       *
       * @param communication_channel Sets an offset to the MPI_Isend and
       * MPI_Irecv calls that avoids interference with other ongoing
       * import_from_ghosted_array_start() calls on different
       * entries. Typically handled within the blocks of a block vector.
       * Any value less than 200 is a valid value.
       *
       * @param ghost_array The array of ghost data that is sent to a remote
       * owner of the respective index in a vector. Its size must either be
       * n_ghost_indices() or equal the number of ghost indices in the larger
       * index set that was given as second argument to
       * set_ghost_indices(). This or the subsequent
       * import_from_ghosted_array_finish() function, the order is
       * implementation-dependent, will set all data entries behind @p
       * ghost_array to zero.
       *
       * @param temporary_storage A temporary storage array of length
       * n_import_indices() that is used to hold the packed data from MPI
       * communication that will later be written into the locally owned
       * array. Note that this array must not be touched until the respective
       * import_from_ghosted_array_finish() call has been made because the
       * model uses non-blocking communication.
       *
       * @param requests The list of MPI requests for the ongoing non-blocking
       * communication that will be finalized in the
       * export_to_ghosted_array_finish() call.
       *
       * This functionality is used in
       * LinearAlgebra::distributed::Vector::compress().
       */
      template <typename Number, typename MemorySpaceType = MemorySpace::Host>
      void
      import_from_ghosted_array_start(
        const VectorOperation::values             vector_operation,
        const unsigned int                        communication_channel,
        const ArrayView<Number, MemorySpaceType> &ghost_array,
        const ArrayView<Number, MemorySpaceType> &temporary_storage,
        std::vector<MPI_Request>                 &requests) const;

      /**
       * Finish importing the data from an array indexed by the ghost
       * indices of this class into a specified locally owned array, combining
       * the results according to the given input @p vector_operation.
       *
       * @param vector_operation Defines how the data sent to the owner should
       * be combined with the existing entries, e.g., added into.
       *
       * @param temporary_storage The same array given to the
       * import_from_ghosted_array_start() call that contains the packed data
       * from MPI communication. In thus function, it is combined at the
       * corresponding entries described by the ghost relations according to
       * @p vector_operation.
       *
       * @param ghost_array The array of ghost data that is sent to a remote
       * owner of the respective index in a vector. Its size must either be
       * n_ghost_indices() or equal the number of ghost indices in the larger
       * index set that was given as second argument to
       * set_ghost_indices(). This function will set all data entries behind
       * @p ghost_array to zero for the implementation-dependent cases when it
       * was not already done in the import_from_ghosted_array_start() call.
       *
       * @param locally_owned_storage The array of data where the resulting data
       * sent by remote processes to the calling process will be accumulated
       * into.
       *
       * @param requests The list of MPI requests for the ongoing non-blocking
       * communication that have been initiated in the
       * import_to_ghosted_array_finish() call. This must be the same array as
       * passed to that function, otherwise MPI will likely throw an error.
       *
       * This functionality is used in
       * LinearAlgebra::distributed::Vector::compress().
       */
      template <typename Number, typename MemorySpaceType = MemorySpace::Host>
      void
      import_from_ghosted_array_finish(
        const VectorOperation::values                   vector_operation,
        const ArrayView<const Number, MemorySpaceType> &temporary_storage,
        const ArrayView<Number, MemorySpaceType>       &locally_owned_storage,
        const ArrayView<Number, MemorySpaceType>       &ghost_array,
        std::vector<MPI_Request>                       &requests) const;
#endif

      /**
       * Compute the memory consumption of this structure.
       */
      std::size_t
      memory_consumption() const;

      /**
       * Exception
       */
      DeclException2(ExcIndexNotPresent,
                     types::global_dof_index,
                     unsigned int,
                     << "Global index " << arg1
                     << " neither owned nor ghost on proc " << arg2 << '.');

      /**
       * Exception
       */
      DeclException3(ExcGhostIndexArrayHasWrongSize,
                     unsigned int,
                     unsigned int,
                     unsigned int,
                     << "The size of the ghost index array (" << arg1
                     << ") must either equal the number of ghost in the "
                     << "partitioner (" << arg2
                     << ") or be equal in size to a more comprehensive index"
                     << "set which contains " << arg3
                     << " elements for this partitioner.");

    private:
      /**
       * Initialize import_indices_plain_dev from import_indices_data. This
       * function is only used when using @ref GlossDevice "device"-aware MPI.
       */
      void
      initialize_import_indices_plain_dev() const;

      /**
       * The global size of the vector over all processors
       */
      types::global_dof_index global_size;

      /**
       * The range of the vector that is stored locally.
       */
      IndexSet locally_owned_range_data;

      /**
       * The range of the vector that is stored locally. Extracted from
       * locally_owned_range for performance reasons.
       */
      std::pair<types::global_dof_index, types::global_dof_index>
        local_range_data;

      /**
       * The set of indices to which we need to have read access but that are
       * not locally owned
       */
      IndexSet ghost_indices_data;

      /**
       * A variable caching the number of ghost indices. It would be expensive
       * to use @p ghost_indices.n_elements() to compute this.
       */
      unsigned int n_ghost_indices_data;

      /**
       * An array that contains information which processors my ghost indices
       * belong to and how many those indices are
       */
      std::vector<std::pair<unsigned int, unsigned int>> ghost_targets_data;

      /**
       * The set of (local) indices that we are importing during compress(),
       * i.e., others' ghosts that belong to the local range. Similar
       * structure as in an IndexSet, but tailored to be iterated over, and
       * some indices may be duplicates.
       */
      std::vector<std::pair<unsigned int, unsigned int>> import_indices_data;

      /**
       * The set of (local) indices that we are importing during compress(),
       * i.e., others' ghosts that belong to the local range. The data stored is
       * the same as in import_indices_data but the data is expanded in plain
       * arrays. This variable is only used when using @ref GlossDevice "device"-aware MPI.
       */
      // The variable is mutable to enable lazy initialization in
      // export_to_ghosted_array_start().
      mutable std::vector<
        Kokkos::View<unsigned int *, MemorySpace::Default::kokkos_space>>
        import_indices_plain_dev;

      /**
       * A variable caching the number of ghost indices. It would be expensive
       * to compute it by iterating over the import indices and accumulate them.
       */
      unsigned int n_import_indices_data;

      /**
       * The set of processors and length of data field which send us their
       * ghost data
       */
      std::vector<std::pair<unsigned int, unsigned int>> import_targets_data;

      /**
       * An array that caches the number of chunks in the import indices per MPI
       * rank. The length is import_indices_data.size()+1.
       */
      std::vector<unsigned int> import_indices_chunks_by_rank_data;

      /**
       * A variable caching the number of ghost indices in a larger set of
       * indices given by the optional argument to set_ghost_indices().
       */
      unsigned int n_ghost_indices_in_larger_set;

      /**
       * An array that caches the number of chunks in the import indices per MPI
       * rank. The length is ghost_indices_subset_data.size()+1.
       */
      std::vector<unsigned int> ghost_indices_subset_chunks_by_rank_data;

      /**
       * The set of indices that appear for an IndexSet that is a subset of a
       * larger set. Similar structure as in an IndexSet within all ghost
       * indices, but tailored to be iterated over.
       */
      std::vector<std::pair<unsigned int, unsigned int>>
        ghost_indices_subset_data;

      /**
       * The ID of the current processor in the MPI network
       */
      unsigned int my_pid;

      /**
       * The total number of processors active in the problem
       */
      unsigned int n_procs;

      /**
       * The MPI communicator involved in the problem
       */
      MPI_Comm communicator;

      /**
       * A variable storing whether the ghost indices have been explicitly set.
       */
      bool have_ghost_indices;
    };



    /*--------------------- Inline functions --------------------------------*/

#ifndef DOXYGEN

    inline types::global_dof_index
    Partitioner::size() const
    {
      return global_size;
    }



    inline const IndexSet &
    Partitioner::locally_owned_range() const
    {
      return locally_owned_range_data;
    }



    inline std::pair<types::global_dof_index, types::global_dof_index>
    Partitioner::local_range() const
    {
      return local_range_data;
    }



    inline unsigned int
    Partitioner::locally_owned_size() const
    {
      types::global_dof_index size =
        local_range_data.second - local_range_data.first;
      Assert(size <= std::numeric_limits<unsigned int>::max(),
             ExcNotImplemented());
      return static_cast<unsigned int>(size);
    }



    inline bool
    Partitioner::in_local_range(
      const types::global_dof_index global_index) const
    {
      return (local_range_data.first <= global_index &&
              global_index < local_range_data.second);
    }



    inline bool
    Partitioner::is_ghost_entry(
      const types::global_dof_index global_index) const
    {
      // if the index is in the global range, it is trivially not a ghost
      if (in_local_range(global_index) == true)
        return false;
      else
        return ghost_indices().is_element(global_index);
    }



    inline unsigned int
    Partitioner::global_to_local(
      const types::global_dof_index global_index) const
    {
      Assert(in_local_range(global_index) || is_ghost_entry(global_index),
             ExcIndexNotPresent(global_index, my_pid));
      if (in_local_range(global_index))
        return static_cast<unsigned int>(global_index - local_range_data.first);
      else
        {
          // avoid checking the ghost index set via a binary search twice by
          // querying the index within set, which returns invalid_dof_index
          // for non-existent entries
          const types::global_dof_index index_within_ghosts =
            ghost_indices_data.index_within_set(global_index);
          if (index_within_ghosts == numbers::invalid_dof_index)
            return numbers::invalid_unsigned_int;
          else
            return locally_owned_size() +
                   static_cast<unsigned int>(index_within_ghosts);
        }
    }



    inline types::global_dof_index
    Partitioner::local_to_global(const unsigned int local_index) const
    {
      AssertIndexRange(local_index,
                       locally_owned_size() + n_ghost_indices_data);
      if (local_index < locally_owned_size())
        return local_range_data.first + types::global_dof_index(local_index);
      else
        return ghost_indices_data.nth_index_in_set(local_index -
                                                   locally_owned_size());
    }



    inline const IndexSet &
    Partitioner::ghost_indices() const
    {
      return ghost_indices_data;
    }



    inline unsigned int
    Partitioner::n_ghost_indices() const
    {
      return n_ghost_indices_data;
    }



    inline const std::vector<std::pair<unsigned int, unsigned int>> &
    Partitioner::ghost_indices_within_larger_ghost_set() const
    {
      return ghost_indices_subset_data;
    }



    inline const std::vector<std::pair<unsigned int, unsigned int>> &
    Partitioner::ghost_targets() const
    {
      return ghost_targets_data;
    }


    inline const std::vector<std::pair<unsigned int, unsigned int>> &
    Partitioner::import_indices() const
    {
      return import_indices_data;
    }



    inline unsigned int
    Partitioner::n_import_indices() const
    {
      return n_import_indices_data;
    }



    inline const std::vector<std::pair<unsigned int, unsigned int>> &
    Partitioner::import_targets() const
    {
      return import_targets_data;
    }



    inline unsigned int
    Partitioner::this_mpi_process() const
    {
      // return the id from the variable stored in this class instead of
      // Utilities::MPI::this_mpi_process() in order to make this query also
      // work when MPI is not initialized.
      return my_pid;
    }



    inline unsigned int
    Partitioner::n_mpi_processes() const
    {
      // return the number of MPI processes from the variable stored in this
      // class instead of Utilities::MPI::n_mpi_processes() in order to make
      // this query also work when MPI is not initialized.
      return n_procs;
    }



    inline MPI_Comm
    Partitioner::get_mpi_communicator() const
    {
      return communicator;
    }



    inline bool
    Partitioner::ghost_indices_initialized() const
    {
      return have_ghost_indices;
    }

#endif // ifndef DOXYGEN

  } // end of namespace MPI

} // end of namespace Utilities


DEAL_II_NAMESPACE_CLOSE

#endif
