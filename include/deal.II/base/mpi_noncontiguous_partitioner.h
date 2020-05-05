// ---------------------------------------------------------------------
//
// Copyright (C) 2020 by the deal.II authors
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

#ifndef dealii_mpi_noncontiguous_partitioner_h
#define dealii_mpi_noncontiguous_partitioner_h

#include <deal.II/base/config.h>

#include <deal.II/base/mpi.h>
#include <deal.II/base/mpi_compute_index_owner_internal.h>
#include <deal.II/base/mpi_tags.h>

#include <deal.II/lac/communication_pattern_base.h>
#include <deal.II/lac/vector_space_vector.h>


DEAL_II_NAMESPACE_OPEN

namespace Utilities
{
  namespace MPI
  {
    /**
     * A flexible Partitioner class, which does not impose restrictions
     * regarding the order of the underlying index sets.
     *
     * @author Peter Munch, 2020
     */
    class NoncontiguousPartitioner
      : public dealii::LinearAlgebra::CommunicationPatternBase
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
                               const MPI_Comm &communicator);

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
        const MPI_Comm &                            communicator);

      /**
       * Fill the vector @p ghost_array according to the precomputed communication
       * pattern with values from @p locally_owned_array.
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
      template <typename Number>
      void
      export_to_ghosted_array(
        const ArrayView<const Number> &locally_owned_array,
        const ArrayView<Number> &      ghost_array) const;

      /**
       * Same as above but with an interface similar to
       * Utilities::MPI::Partitioner::export_to_ghosted_array_start and
       * Utilities::MPI::Partitioner::export_to_ghosted_array_finish. In this
       * function, the user can provide the temporary data structures to be
       * used.
       *
       * @pre The size of the @p temporary_storage vector has to be at least
       *   as large as the sum of the number of entries in the index sets
       *   passed to the constructor and the reinit() functions. The reason
       *   for this is that this vector is used as buffer for both sending
       *   and receiving data.
       */
      template <typename Number>
      void
      export_to_ghosted_array(
        const unsigned int             communication_channel,
        const ArrayView<const Number> &locally_owned_array,
        const ArrayView<Number> &      temporary_storage,
        const ArrayView<Number> &      ghost_array,
        std::vector<MPI_Request> &     requests) const;

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
       * above.
       */
      template <typename Number>
      void
      export_to_ghosted_array_start(
        const unsigned int             communication_channel,
        const ArrayView<const Number> &locally_owned_array,
        const ArrayView<Number> &      temporary_storage,
        std::vector<MPI_Request> &     requests) const;

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
       * above.
       */
      template <typename Number>
      void
      export_to_ghosted_array_finish(
        const ArrayView<const Number> &temporary_storage,
        const ArrayView<Number> &      ghost_array,
        std::vector<MPI_Request> &     requests) const;

      /**
       * Returns the number of processes this process sends data to and the
       * number of processes this process receives data from.
       */
      std::pair<unsigned int, unsigned int>
      n_targets();

      /**
       * Return memory consumption in Byte.
       */
      types::global_dof_index
      memory_consumption();

      /**
       * Return the underlying communicator.
       */
      const MPI_Comm &
      get_mpi_communicator() const override;

      /**
       * Initialize the inner data structures.
       */
      void
      reinit(const IndexSet &indexset_locally_owned,
             const IndexSet &indexset_ghost,
             const MPI_Comm &communicator) override;

      /**
       * Initialize the inner data structures.
       */
      void
      reinit(const std::vector<types::global_dof_index> &indices_locally_owned,
             const std::vector<types::global_dof_index> &indices_ghost,
             const MPI_Comm &                            communicator);

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
       * @note Together with `recv_ptr` this forms a CRS data structure.
       */
      std::vector<types::global_dof_index> recv_indices;

      /**
       * Buffer containing the values sorted by rank for sending and receiving.
       *
       * @note Only allocated if not provided externally by user.
       *
       * @note At this place we do not know the type of the data to be sent. So
       *   we use an arbitrary type of size 1 byte. The type is cast to the
       *   requested type in the relevant functions.
       */
      mutable std::vector<uint8_t> buffers;

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
