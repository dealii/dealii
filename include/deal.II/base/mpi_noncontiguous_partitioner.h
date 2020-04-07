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

#ifndef dealii_mpi_noncontiguous_vector_h
#define dealii_mpi_noncontiguous_vector_h

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
    template <typename Number = double>
    class NoncontiguousPartitioner
      : public LinearAlgebra::CommunicationPatternBase
    {
    public:
      /**
       * Default constructor. Requires calling one of the reinit() functions
       * to create a valid object.
       */
      NoncontiguousPartitioner() = default;

      /**
       * Constructor. Set up point-to-point communication pattern based on the
       * IndexSets arguments @p indexset_has and @p indexset_want for the MPI
       * communicator @p communicator.
       */
      NoncontiguousPartitioner(const IndexSet &indexset_has,
                               const IndexSet &indexset_want,
                               const MPI_Comm &communicator);

      /**
       * Constructor. Same as above but for vectors of indices @p indices_has
       * and @p indices_want. This allows the indices to not be sorted and the
       * values are read and written automatically at the right position of
       * the vector during update_values(), update_values_start(), and
       * update_values_finish(). It is allowed to include entries with the
       * value numbers::invalid_dof_index which do not take part of the index
       * exchange but are present in the data vectors as padding.
       */
      NoncontiguousPartitioner(
        const std::vector<types::global_dof_index> &indices_has,
        const std::vector<types::global_dof_index> &indices_want,
        const MPI_Comm &                            communicator);

      /**
       * Fill the vector @p dst according to the precomputed communication
       * pattern with values from @p src.
       *
       * @pre The vectors only have to provide a method begin(), which allows
       *   to access their raw data.
       *
       * @note This function calls the methods update_values_start() and
       *   update_values_finish() in sequence. Users can call these two
       *   functions separately and hereby overlap communication and
       *   computation.
       */
      template <typename VectorType>
      void
      update_values(VectorType &dst, const VectorType &src) const;

      /**
       * Start update: Data is packed, non-blocking send and receives
       * are started.
       */
      template <typename VectorType>
      void
      update_values_start(const VectorType &src, const unsigned int tag) const;

      /**
       * Finish update. The method waits until all data has been sent and
       * received. Once data from any process is received it is processed and
       * placed at the right position of the vector @p dst.
       */
      template <typename VectorType>
      void
      update_values_finish(VectorType &dst, const unsigned int tag) const;

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
      reinit(const IndexSet &indexset_has,
             const IndexSet &indexset_want,
             const MPI_Comm &communicator) override;

      /**
       * Initialize the inner data structures.
       */
      void
      reinit(const std::vector<types::global_dof_index> &indices_has,
             const std::vector<types::global_dof_index> &indices_want,
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
       *  Buffer containing the values sorted accoding to the ranks.
       */
      mutable std::vector<Number> send_buffers;

      /**
       * MPI requests for sending.
       */
      mutable std::vector<MPI_Request> send_requests;

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
       * Buffer containing the values sorted by rank.
       */
      mutable std::vector<Number> recv_buffers;

      /**
       * MPI requests for receiving.
       */
      mutable std::vector<MPI_Request> recv_requests;
    };

  } // namespace MPI
} // namespace Utilities

DEAL_II_NAMESPACE_CLOSE

#endif
