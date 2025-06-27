// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_mpi_noncontiguous_partitioner_h
#define dealii_mpi_noncontiguous_partitioner_h

#include <deal.II/base/config.h>

#include <deal.II/base/array_view.h>
#include <deal.II/base/communication_pattern_base.h>
#include <deal.II/base/mpi_stub.h>
#include <deal.II/base/types.h>

#include <deal.II/lac/vector_operation.h>


DEAL_II_NAMESPACE_OPEN

namespace Utilities
{
  namespace MPI
  {
    /**
     * A flexible Partitioner class, which does not impose restrictions
     * regarding the order of the underlying index sets. In other words,
     * this class implements the interface of the
     * Utilities::MPI::CommunicationPatternBase base class with no
     * assumption that every process stores a contiguous part of the
     * array of objects, but that indeed the locally owned indices
     * can be an arbitrary subset of all indices of elements of the array
     * to which they refer.
     *
     * If you want to store only contiguous parts of these arrays on
     * each process, take a look at Utilities::MPI::Partitioner.
     */
    class NoncontiguousPartitioner
      : public Utilities::MPI::CommunicationPatternBase
    {
    public:
      /**
       * Default constructor. Requires calling one of the reinit() functions
       * to create a valid object.
       */
      NoncontiguousPartitioner() = default;

      /**
       * Constructor. Set up point-to-point communication pattern based on the
       * IndexSets arguments @p indexset_locally_owned and @p indexset_ghost for the MPI
       * communicator @p communicator.
       */
      NoncontiguousPartitioner(const IndexSet &indexset_locally_owned,
                               const IndexSet &indexset_ghost,
                               const MPI_Comm  communicator);

      /**
       * Constructor. Same as above but for vectors of indices @p indices_locally_owned
       * and @p indices_ghost. This allows the indices to not be sorted and the
       * values are read and written automatically at the right position of
       * the vector during update_values(), update_values_start(), and
       * update_values_finish(). It is allowed to include entries with the
       * value numbers::invalid_dof_index which do not take part of the index
       * exchange but are present in the data vectors as padding.
       */
      NoncontiguousPartitioner(
        const std::vector<types::global_dof_index> &indices_locally_owned,
        const std::vector<types::global_dof_index> &indices_ghost,
        const MPI_Comm                              communicator);

      /**
       * Fill the vector @p ghost_array according to the precomputed communication
       * pattern with values from @p locally_owned_array.
       *
       * In the default case, only one object is communicated per entry
       * (`n_components_templated == 1'). If you want to communicate more
       * entries, you can increase the value of @p n_components_templated in the
       * case that you know the size at compile time. If you want to set the
       * size during runtime, you can set @p n_components. However,
       * @p n_components_templated has to be set to `0` in this case. Either
       * @p n_components_templated or @p n_components can be set.
       *
       * @pre The vectors only have to provide a method begin(), which allows
       *   to access their raw data.
       *
       * @pre The size of both vectors must be at least as large as the number
       *   of entries in the index sets passed to the constructors or the
       *   reinit() functions.
       *
       * @note This function calls the methods update_values_start() and
       *   update_values_finish() in sequence. Users can call these two
       *   functions separately and hereby overlap communication and
       *   computation.
       */
      template <typename Number, unsigned int n_components_templated = 1>
      void
      export_to_ghosted_array(
        const ArrayView<const Number> &locally_owned_array,
        const ArrayView<Number>       &ghost_array,
        const unsigned int             n_components = 0) const;

      /**
       * Same as above but with an interface similar to
       * Utilities::MPI::Partitioner::export_to_ghosted_array_start and
       * Utilities::MPI::Partitioner::export_to_ghosted_array_finish. In this
       * function, the user can provide the temporary data structures to be
       * used.
       *
       * @pre The size of the @p temporary_storage vector has to be at least
       *   temporary_storage_size. The reason for this is that this vector is
       *   used as buffer for both sending and receiving data.
       *
       * @note Any value less than 10 is a valid value of
       *   @p communication_channel.
       */
      template <typename Number, unsigned int n_components_templated = 1>
      void
      export_to_ghosted_array(
        const unsigned int             communication_channel,
        const ArrayView<const Number> &locally_owned_array,
        const ArrayView<Number>       &temporary_storage,
        const ArrayView<Number>       &ghost_array,
        std::vector<MPI_Request>      &requests,
        const unsigned int             n_components = 0) const;

      /**
       * Start update: Data is packed, non-blocking send and receives
       * are started.
       *
       * @note In contrast to the function
       *   Utilities::MPI::Partitioner::export_to_ghosted_array_start, the user
       *   does not pass a reference to the destination vector, since the data
       *   is received into a designated part of the buffer @p temporary_storage. This
       *   allows for padding and other post-processing of the received data.
       *
       * @pre The required size of the vectors are the same as in the functions
       *   above.
       *
       * @note Any value less than 10 is a valid value of
       *   @p communication_channel.
       */
      template <typename Number, unsigned int n_components_templated = 1>
      void
      export_to_ghosted_array_start(
        const unsigned int             communication_channel,
        const ArrayView<const Number> &locally_owned_array,
        const ArrayView<Number>       &temporary_storage,
        std::vector<MPI_Request>      &requests,
        const unsigned int             n_components = 0) const;

      /**
       * Finish update. The method waits until all data has been sent and
       * received. Once data from any process is received it is processed and
       * placed at the right position of the vector @p dst.
       *
       * @note In contrast to the function
       *   Utilities::MPI::Partitioner::export_to_ghosted_array_finish, the user
       *   also has to pass a reference to the buffer @p temporary_storage,
       *   since the data has been received into the buffer and not into the
       *   destination vector.
       *
       * @pre The required size of the vectors are the same as in the functions
       *   above.
       */
      template <typename Number, unsigned int n_components_templated = 1>
      void
      export_to_ghosted_array_finish(
        const ArrayView<const Number> &temporary_storage,
        const ArrayView<Number>       &ghost_array,
        std::vector<MPI_Request>      &requests,
        const unsigned int             n_components = 0) const;

      /**
       * Similar to the above functions but for importing vector entries
       * from @p ghost_array to @p locally_owned_storage.
       *
       * @note In contrast to the functions in
       *   Utilities::MPI::Partitioner, this function expects that
       *   locally_owned_storage is empty.
       */
      template <typename Number>
      void
      import_from_ghosted_array(
        const VectorOperation::values vector_operation,
        const ArrayView<Number>      &ghost_array,
        const ArrayView<Number>      &locally_owned_storage) const;

      /**
       * Similar to the above function with the difference that
       * users can provide temporaty arrays. This function calls
       * import_from_ghosted_array_start() and
       * import_from_ghosted_array_finish() in sequence.
       */
      template <typename Number>
      void
      import_from_ghosted_array(const VectorOperation::values vector_operation,
                                const unsigned int        communication_channel,
                                const ArrayView<Number>  &ghost_array,
                                const ArrayView<Number>  &temporary_storage,
                                const ArrayView<Number>  &locally_owned_storage,
                                std::vector<MPI_Request> &requests) const;

      /**
       * Start update for importig values: Data is packed, non-blocking send
       * and receives are started.
       */
      template <typename Number>
      void
      import_from_ghosted_array_start(
        const VectorOperation::values vector_operation,
        const unsigned int            communication_channel,
        const ArrayView<Number>      &ghost_array,
        const ArrayView<Number>      &temporary_storage,
        std::vector<MPI_Request>     &requests) const;

      /**
       * Finish update for importing values. The method waits until all data has
       * been sent and received. Once data from any process is received it is
       * processed and placed at the right position of the vector
       * @p locally_owned_storage.
       */
      template <typename Number>
      void
      import_from_ghosted_array_finish(
        const VectorOperation::values  vector_operation,
        const ArrayView<const Number> &temporary_storage,
        const ArrayView<Number>       &locally_owned_storage,
        std::vector<MPI_Request>      &requests) const;

      /**
       * Returns the number of processes this process sends data to and the
       * number of processes this process receives data from.
       */
      std::pair<unsigned int, unsigned int>
      n_targets() const;

      /**
       * Return the size of the temporary storage needed by the
       * export_to_ghosted_array() functions, if the temporary storage is
       * handled by the user code.
       */
      unsigned int
      temporary_storage_size() const;

      /**
       * Return memory consumption in Byte.
       */
      types::global_dof_index
      memory_consumption();

      /**
       * Return the underlying MPI communicator.
       */
      MPI_Comm
      get_mpi_communicator() const override;

      void
      reinit(const IndexSet &locally_owned_indices,
             const IndexSet &ghost_indices,
             const MPI_Comm  communicator) override;

      /**
       * Initialize the inner data structures using explicit sets of
       * indices. See the documentation of the other reinit() function for
       * what the function does.
       */
      void
      reinit(const std::vector<types::global_dof_index> &locally_owned_indices,
             const std::vector<types::global_dof_index> &ghost_indices,
             const MPI_Comm                              communicator);

    private:
      /**
       * MPI communicator.
       */
      MPI_Comm communicator;

      /**
       * The ranks this process sends data to.
       */
      std::vector<unsigned int> send_ranks;

      /**
       * Offset of each process within send_buffer.
       *
       * @note Together with `send_indices` this forms a CRS data structure.
       */
      std::vector<types::global_dof_index> send_ptr;

      /**
       * Local index of each entry in send_buffer within the destination
       * vector.
       *
       * @note Together with `send_ptr` this forms a CRS data structure.
       */
      std::vector<types::global_dof_index> send_indices;

      /**
       * The ranks this process receives data from.
       */
      std::vector<unsigned int> recv_ranks;

      /**
       * Offset of each process within recv_buffer.
       *
       * @note Together with `recv_indices` this forms a CRS data structure.
       */
      std::vector<types::global_dof_index> recv_ptr;

      /**
       * Local index of each entry in recv_buffer within the destination
       * vector.
       *
       * In the case that an entry is more than once requested by the
       * same rank, local index of each entry in recv_buffer within
       * recv_indices_duplicates_ptr.
       *
       * @note Together with `recv_ptr` this forms a CRS data structure.
       */
      std::vector<types::global_dof_index> recv_indices;

      /**
       * In the case that an entry is more than once requested by the
       * same rank, offset of each index in recv_indices_duplicates.
       */
      std::vector<unsigned int> recv_indices_duplicates_ptr;

      /**
       * Local index of each entry within the destination vector in the
       * case that an entry is more than once requested by the same rank.
       *
       * @note Together with `recv_indices_duplicates_ptr`
       * this forms a CRS data structure.
       */
      std::vector<unsigned int> recv_indices_duplicates;

      /**
       * Buffer containing the values sorted by rank for sending and receiving.
       *
       * @note Only allocated if not provided externally by user.
       *
       * @note At this place we do not know the type of the data to be sent. So
       *   we use an arbitrary type of size 1 byte. The type is cast to the
       *   requested type in the relevant functions.
       */
      mutable std::vector<std::uint8_t> buffers;

      /**
       * MPI requests for sending and receiving.
       *
       * @note Only allocated if not provided externally by user.
       */
      mutable std::vector<MPI_Request> requests;
    };

  } // namespace MPI
} // namespace Utilities

DEAL_II_NAMESPACE_CLOSE

#endif
