// ---------------------------------------------------------------------
//
// Copyright (C) 2019 - 2022 by the deal.II authors
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

#ifndef dealii_base_mpi_compute_index_owner_internal_h
#define dealii_base_mpi_compute_index_owner_internal_h

#include <deal.II/base/config.h>

#include <deal.II/base/mpi.h>
#include <deal.II/base/mpi_consensus_algorithms.h>

DEAL_II_NAMESPACE_OPEN

namespace Utilities
{
  namespace MPI
  {
    namespace internal
    {
      /**
       * An internal namespace used for Utilities::MPI::compute_index_owner()
       * and for Utilities::MPI::Partitioner::set_ghost_indices().
       */
      namespace ComputeIndexOwner
      {
        class FlexibleIndexStorage
        {
        public:
          using index_type = unsigned int;
          static const index_type invalid_index_value =
            numbers::invalid_unsigned_int;

          FlexibleIndexStorage(const bool use_vector = true);

          void
          reinit(const bool        use_vector,
                 const bool        index_range_contiguous,
                 const std::size_t size);

          void
          fill(const std::size_t start,
               const std::size_t end,
               const index_type &value);

          index_type &
          operator[](const std::size_t index);

          index_type
          operator[](const std::size_t index) const;

          bool
          entry_has_been_set(const std::size_t index) const;

        private:
          bool                              use_vector;
          std::size_t                       size;
          std::vector<index_type>           data;
          std::map<std::size_t, index_type> data_map;
        };

        /**
         * Specialization of ConsensusAlgorithms::Process for setting up the
         * Dictionary even if there are ranges in the IndexSet space not owned
         * by any processes.
         *
         * @note Only for internal usage.
         */
        class DictionaryPayLoad
          : public ConsensusAlgorithms::Process<
              std::vector<
                std::pair<types::global_dof_index, types::global_dof_index>>,
              std::vector<unsigned int>>
        {
        public:
          /**
           * Constructor.
           */
          DictionaryPayLoad(
            const std::map<unsigned int,
                           std::vector<std::pair<types::global_dof_index,
                                                 types::global_dof_index>>>
              &                   buffers,
            FlexibleIndexStorage &actually_owning_ranks,
            const std::pair<types::global_dof_index, types::global_dof_index>
              &                        local_range,
            std::vector<unsigned int> &actually_owning_rank_list);

          /**
           * Implementation of
           * Utilities::MPI::ConsensusAlgorithms::Process::compute_targets().
           */
          virtual std::vector<unsigned int>
          compute_targets() override;

          /**
           * Implementation of
           * Utilities::MPI::ConsensusAlgorithms::Process::create_request().
           */
          virtual void
          create_request(const unsigned int other_rank,
                         std::vector<std::pair<types::global_dof_index,
                                               types::global_dof_index>>
                           &send_buffer) override;

          /**
           * Implementation of
           * Utilities::MPI::ConsensusAlgorithms::Process::answer_request().
           */
          virtual void
          answer_request(
            const unsigned int                                     other_rank,
            const std::vector<std::pair<types::global_dof_index,
                                        types::global_dof_index>> &buffer_recv,
            std::vector<unsigned int> &request_buffer) override;

        private:
          const std::map<unsigned int,
                         std::vector<std::pair<types::global_dof_index,
                                               types::global_dof_index>>>
            &buffers;

          FlexibleIndexStorage &actually_owning_ranks;

          const std::pair<types::global_dof_index, types::global_dof_index>
            &local_range;

          std::vector<unsigned int> &actually_owning_rank_list;
        };



        /**
         * Dictionary class with basic partitioning in terms of a single
         * interval of fixed size known to all MPI ranks for two-stage index
         * lookup.
         */
        struct Dictionary
        {
          /**
           * The minimum grain size for the intervals.
           *
           * We choose to limit the smallest size an interval for the
           * two-stage lookup can have with the following two conflicting
           * goals in mind: On the one hand, we do not want intervals in the
           * dictionary to become too short. For uneven distributions of
           * unknowns (some ranks with several thousands of unknowns, others
           * with none), the lookup DoFs -> dictionary then involves sending
           * from one MPI rank to many other MPI ranks holding dictionary
           * intervals, leading to an exceedingly high number of messages some
           * ranks have to send. Also, fewer longer intervals are generally
           * more efficient to look up. On the other hand, a range size too
           * large leads to opposite effect of many messages that come into a
           * particular dictionary owner in the lookup DoFs ->
           * dictionary. With the current setting, we get at most 64 messages
           * coming to a single MPI rank in case there is 1 dof per MPI rank,
           * which is reasonably low. At the same time, uneven distributions
           * up to factors of 4096 can be handled with at most 64 messages as
           * well.
           */
          static constexpr unsigned int range_minimum_grain_size = 64;

          /**
           * Factor that determines if an index set is sparse or not. An index
           * set if sparse if less than 25% of the indices are owned by any
           * process. If the index set is sparse, we switch the internal storage
           * from a fast storage (vector) to a memory-efficient storage (map).
           */
          static constexpr unsigned int sparsity_factor = 4;

          /**
           * A vector with as many entries as there are dofs in the dictionary
           * of the current process, and each entry containing the rank of the
           * owner of that dof in the IndexSet `owned_indices`. This is
           * queried in the index lookup, so we keep an expanded list.
           */
          FlexibleIndexStorage actually_owning_ranks;

          /**
           * A sorted vector containing the MPI ranks appearing in
           * `actually_owning_ranks`.
           */
          std::vector<unsigned int> actually_owning_rank_list;

          /**
           * The number of unknowns in the dictionary for on each MPI rank
           * used for the index space splitting. For simplicity of index
           * lookup without additional communication, this number is the same
           * on all MPI ranks.
           */
          types::global_dof_index dofs_per_process;

          /**
           * The local range of the global index space that is represented in
           * the dictionary, computed from `dofs_per_process`, the current
           * MPI rank, and range_minimum_grain_size.
           */
          std::pair<types::global_dof_index, types::global_dof_index>
            local_range;

          /**
           * The actual size, computed as the minimum of dofs_per_process and
           * the possible end of the index space. Equivalent to
           * `local_range.second - local_range.first`.
           */
          types::global_dof_index locally_owned_size;

          /**
           * The global size of the index space.
           */
          types::global_dof_index size;

          /**
           * The number of ranks the `owned_indices` IndexSet is distributed
           * among.
           */
          unsigned int n_dict_procs_in_owned_indices;

          /**
           * A stride to distribute the work more evenly over MPI ranks in
           * case the grain size forces us to have fewer ranges than we have
           * processors.
           */
          unsigned int stride_small_size;

          /**
           * Set up the dictionary by computing the partitioning from the
           * global size and sending the rank information on locally owned
           * ranges to the owner of the dictionary part.
           */
          void
          reinit(const IndexSet &owned_indices, const MPI_Comm &comm);

          /**
           * Translate a global dof index to the MPI rank in the dictionary
           * using `dofs_per_process`. We multiply by `stride_small_size` to
           * ensure a balance over the MPI ranks due to the grain size.
           */
          unsigned int
          dof_to_dict_rank(const types::global_dof_index i);

          /**
           * Given an MPI rank id of an arbitrary processor, return the index
           * offset where the local range of that processor begins.
           */
          types::global_dof_index
          get_index_offset(const unsigned int rank);

          /**
           * Given the rank in the owned indices from `actually_owning_ranks`,
           * this returns the index of the rank in the
           * `actually_owning_rank_list`.
           */
          unsigned int
          get_owning_rank_index(const unsigned int rank_in_owned_indices,
                                const unsigned int guess = 0);

        private:
          /**
           * Compute the partition from the global size of the index space and
           * the number of ranks.
           */
          void
          partition(const IndexSet &owned_indices, const MPI_Comm &comm);
        };



        /**
         * Specialization of ConsensusAlgorithms::Process for the context of
         * Utilities::MPI::compute_index_owner() and
         * Utilities::MPI::Partitioner::set_ghost_indices() with additional
         * payload.
         */
        class ConsensusAlgorithmsPayload
          : public ConsensusAlgorithms::Process<
              std::vector<
                std::pair<types::global_dof_index, types::global_dof_index>>,
              std::vector<unsigned int>>
        {
        public:
          /**
           * Constructor.
           */
          ConsensusAlgorithmsPayload(const IndexSet &owned_indices,
                                     const IndexSet &indices_to_look_up,
                                     const MPI_Comm &comm,
                                     std::vector<unsigned int> &owning_ranks,
                                     const bool track_index_requests = false);

          /**
           * The index space which describes the locally owned space.
           */
          const IndexSet &owned_indices;

          /**
           * The indices which are "ghosts" on a given rank and should be
           * looked up in terms of their owner rank from owned_indices.
           */
          const IndexSet &indices_to_look_up;

          /**
           * The underlying MPI communicator.
           */
          const MPI_Comm comm;

          /**
           * The present MPI rank.
           */
          const unsigned int my_rank;

          /**
           * The total number of ranks participating in the MPI communicator
           * `comm`.
           */
          const unsigned int n_procs;

          /**
           * Controls whether the origin of ghost owner should also be
           * stored. If true, it will be added into `requesters` and can be
           * queried by `get_requesters()`.
           */
          const bool track_index_requests;

          /**
           * The result of the index owner computation: To each index
           * contained in `indices_to_look_up`, this vector contains the MPI
           * rank of the owner in `owned_indices`.
           */
          std::vector<unsigned int> &owning_ranks;

          /**
           * Keeps track of the origin of the requests. The layout of the data
           * structure is as follows: The outermost vector has as many entries
           * as Dictionary::actually_owning_rank_list and represents the
           * information we should send back to the owners from the present
           * dictionary entry. The second vector then collects a list of MPI
           * ranks that have requested data, using the rank in the first pair
           * entry and a list of index ranges as the second entry.
           */
          std::vector<std::vector<
            std::pair<unsigned int,
                      std::vector<std::pair<unsigned int, unsigned int>>>>>
            requesters;

          /**
           * The dictionary handling the requests.
           */
          Dictionary dict;

          /**
           * Array to collect the indices to look up, sorted by the rank in
           * the dictionary.
           */
          std::map<unsigned int, std::vector<types::global_dof_index>>
            indices_to_look_up_by_dict_rank;

          /**
           * The field where the indices for incoming data from the process
           * are stored.
           */
          std::map<unsigned int, std::vector<unsigned int>> recv_indices;

          /**
           * Implementation of
           * Utilities::MPI::ConsensusAlgorithms::Process::answer_request(),
           * adding the owner of a particular index in request_buffer (and
           * keeping track of who requested a particular index in case that
           * information is also desired).
           */
          virtual void
          answer_request(
            const unsigned int                                     other_rank,
            const std::vector<std::pair<types::global_dof_index,
                                        types::global_dof_index>> &buffer_recv,
            std::vector<unsigned int> &request_buffer) override;

          /**
           * Implementation of
           * Utilities::MPI::ConsensusAlgorithms::Process::compute_targets().
           */
          virtual std::vector<unsigned int>
          compute_targets() override;

          /**
           * Implementation of
           * Utilities::MPI::ConsensusAlgorithms::Process::create_request().
           */
          virtual void
          create_request(const unsigned int other_rank,
                         std::vector<std::pair<types::global_dof_index,
                                               types::global_dof_index>>
                           &send_buffer) override;

          /**
           * Implementation of
           * Utilities::MPI::ConsensusAlgorithms::Process::read_answer().
           */
          virtual void
          read_answer(const unsigned int               other_rank,
                      const std::vector<unsigned int> &recv_buffer) override;

          /**
           * Resolve the origin of the requests by sending the information
           * accumulated in terms of the dictionary owners during the run of
           * the consensus algorithm back to the owner in the original
           * IndexSet. This requires some point-to-point communication.
           *
           * @return Map of processors and associated ranges of indices that
           *         are requested from the current rank
           */
          std::map<unsigned int, IndexSet>
          get_requesters();

        private:
          /**
           * Stores the index request in the `requesters` field. We first find
           * out the owner of the index that was requested (using the guess in
           * `owner_index`, as we typically might look up on the same rank
           * several times in a row, which avoids the binary search in
           * Dictionary::get_owning_rank_index()). Once we know the rank of
           * the owner, we fill the vector entry with the rank of the
           * request. Here, we utilize the fact that requests are processed
           * rank-by-rank, so we can simply look at the end of the vector
           * whether there is already some data stored or not. Finally, we
           * build ranges, again using that the index list is sorted and we
           * therefore only need to append at the end.
           */
          void
          append_index_origin(const types::global_dof_index index,
                              unsigned int &                owner_index,
                              const unsigned int            rank_of_request);
        };

        /* ------------------------- inline functions ----------------------- */

        inline unsigned int
        Dictionary::dof_to_dict_rank(const types::global_dof_index i)
        {
          // note: this formula is also explicitly used in
          // get_index_offset(), so keep the two in sync
          return (i / dofs_per_process) * stride_small_size;
        }


        inline types::global_dof_index
        Dictionary::get_index_offset(const unsigned int rank)
        {
          return std::min(dofs_per_process *
                            static_cast<types::global_dof_index>(
                              (rank + stride_small_size - 1) /
                              stride_small_size),
                          size);
        }



        inline unsigned int
        Dictionary::get_owning_rank_index(
          const unsigned int rank_in_owned_indices,
          const unsigned int guess)
        {
          AssertIndexRange(guess, actually_owning_rank_list.size());
          if (actually_owning_rank_list[guess] == rank_in_owned_indices)
            return guess;
          else
            {
              auto it = std::lower_bound(actually_owning_rank_list.begin(),
                                         actually_owning_rank_list.end(),
                                         rank_in_owned_indices);
              Assert(it != actually_owning_rank_list.end(), ExcInternalError());
              Assert(*it == rank_in_owned_indices, ExcInternalError());
              return it - actually_owning_rank_list.begin();
            }
        }


      } // namespace ComputeIndexOwner
    }   // namespace internal
  }     // namespace MPI
} // namespace Utilities

DEAL_II_NAMESPACE_CLOSE

#endif
