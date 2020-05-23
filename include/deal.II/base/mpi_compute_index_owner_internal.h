// ---------------------------------------------------------------------
//
// Copyright (C) 2019 - 2020 by the deal.II authors
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
        /**
         * Specialization of ConsensusAlgorithms::Process for setting up the
         * Dictionary even if there are ranges in the IndexSet space not owned
         * by any processes.
         *
         * @note Only for internal usage.
         */
        class DictionaryPayLoad
          : public ConsensusAlgorithms::Process<
              std::pair<types::global_dof_index, types::global_dof_index>,
              unsigned int>
        {
        public:
          /**
           * Constructor.
           */
          DictionaryPayLoad(
            const std::map<unsigned int,
                           std::vector<std::pair<types::global_dof_index,
                                                 types::global_dof_index>>>
              &                        buffers,
            std::vector<unsigned int> &actually_owning_ranks,
            const std::pair<types::global_dof_index, types::global_dof_index>
              &                        local_range,
            std::vector<unsigned int> &actually_owning_rank_list)
            : buffers(buffers)
            , actually_owning_ranks(actually_owning_ranks)
            , local_range(local_range)
            , actually_owning_rank_list(actually_owning_rank_list)
          {}

          /**
           * Implementation of
           * Utilities::MPI::ConsensusAlgorithms::Process::compute_targets().
           */
          virtual std::vector<unsigned int>
          compute_targets() override
          {
            std::vector<unsigned int> targets;
            for (const auto &rank_pair : buffers)
              targets.push_back(rank_pair.first);

            return targets;
          }

          /**
           * Implementation of
           * Utilities::MPI::ConsensusAlgorithms::Process::create_request().
           */
          virtual void
          create_request(const unsigned int other_rank,
                         std::vector<std::pair<types::global_dof_index,
                                               types::global_dof_index>>
                           &send_buffer) override
          {
            send_buffer = this->buffers.at(other_rank);
          }

          /**
           * Implementation of
           * Utilities::MPI::ConsensusAlgorithms::Process::answer_request().
           */
          virtual void
          answer_request(
            const unsigned int                                     other_rank,
            const std::vector<std::pair<types::global_dof_index,
                                        types::global_dof_index>> &buffer_recv,
            std::vector<unsigned int> &request_buffer) override
          {
            (void)request_buffer; // not needed


            // process message: loop over all intervals
            for (auto interval : buffer_recv)
              {
#ifdef DEBUG
                for (types::global_dof_index i = interval.first;
                     i < interval.second;
                     i++)
                  Assert(actually_owning_ranks[i - local_range.first] ==
                           numbers::invalid_unsigned_int,
                         ExcInternalError());
                Assert(interval.first >= local_range.first &&
                         interval.first < local_range.second,
                       ExcInternalError());
                Assert(interval.second > local_range.first &&
                         interval.second <= local_range.second,
                       ExcInternalError());
#endif
                std::fill(actually_owning_ranks.data() + interval.first -
                            local_range.first,
                          actually_owning_ranks.data() + interval.second -
                            local_range.first,
                          other_rank);
              }
            actually_owning_rank_list.push_back(other_rank);
          }

        private:
          const std::map<unsigned int,
                         std::vector<std::pair<types::global_dof_index,
                                               types::global_dof_index>>>
            &buffers;

          std::vector<unsigned int> &actually_owning_ranks;

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
           * The minimum grain size for the ranges.
           */
          static const unsigned int range_minimum_grain_size = 4096;

          /**
           * A vector with as many entries as there are dofs in the dictionary
           * of the current process, and each entry containing the rank of the
           * owner of that dof in the IndexSet `owned_indices`. This is
           * queried in the index lookup, so we keep an expanded list.
           */
          std::vector<unsigned int> actually_owning_ranks;

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
          types::global_dof_index local_size;

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
          reinit(const IndexSet &owned_indices, const MPI_Comm &comm)
          {
            // 1) set up the partition
            this->partition(owned_indices, comm);

#ifdef DEAL_II_WITH_MPI
            unsigned int my_rank = this_mpi_process(comm);

            types::global_dof_index dic_local_received = 0;
            std::map<unsigned int,
                     std::vector<std::pair<types::global_dof_index,
                                           types::global_dof_index>>>
              buffers;

            std::fill(actually_owning_ranks.begin(),
                      actually_owning_ranks.end(),
                      numbers::invalid_subdomain_id);

            // 2) collect relevant processes and process local dict entries
            for (auto interval = owned_indices.begin_intervals();
                 interval != owned_indices.end_intervals();
                 ++interval)
              {
                // Due to the granularity of the dictionary, the interval
                // might be split into several ranges of processor owner
                // ranks. Here, we process the interval by breaking into
                // smaller pieces in terms of the dictionary number.
                std::pair<types::global_dof_index, types::global_dof_index>
                                   index_range(*interval->begin(), interval->last() + 1);
                const unsigned int owner_last =
                  dof_to_dict_rank(interval->last());
                unsigned int owner_first = numbers::invalid_unsigned_int;
                while (owner_first != owner_last)
                  {
                    Assert(index_range.first < index_range.second,
                           ExcInternalError());

                    owner_first = dof_to_dict_rank(index_range.first);

                    // this explicitly picks up the formula of
                    // dof_to_dict_rank, so the two places must be in sync
                    types::global_dof_index next_index =
                      std::min(get_index_offset(owner_first + 1),
                               index_range.second);

                    Assert(next_index > index_range.first, ExcInternalError());

#  ifdef DEBUG
                    // make sure that the owner is the same on the current
                    // interval
                    for (types::global_dof_index i = index_range.first + 1;
                         i < next_index;
                         ++i)
                      AssertDimension(owner_first, dof_to_dict_rank(i));
#  endif

                    // add the interval, either to the local range or into a
                    // buffer to be sent to another processor
                    if (owner_first == my_rank)
                      {
                        std::fill(actually_owning_ranks.data() +
                                    index_range.first - local_range.first,
                                  actually_owning_ranks.data() + next_index -
                                    local_range.first,
                                  my_rank);
                        dic_local_received += next_index - index_range.first;
                        if (actually_owning_rank_list.empty())
                          actually_owning_rank_list.push_back(my_rank);
                      }
                    else
                      buffers[owner_first].emplace_back(index_range.first,
                                                        next_index);

                    index_range.first = next_index;
                  }
              }

            n_dict_procs_in_owned_indices = buffers.size();
            std::vector<MPI_Request> request;

            // Check if index set space is partitioned globally without gaps.
            if (Utilities::MPI::sum(owned_indices.n_elements(), comm) ==
                owned_indices.size())
              {
                // no gaps: setup is simple! Processes send their locally owned
                // indices to the dictionary. The dictionary stores the sending
                // rank for each index. The dictionary knows exactly
                // when it is set up when all indices it is responsible for
                // have been processed.

                request.reserve(n_dict_procs_in_owned_indices);

                // protect the following communication steps using a mutex:
                static CollectiveMutex      mutex;
                CollectiveMutex::ScopedLock lock(mutex, comm);

                const int mpi_tag =
                  Utilities::MPI::internal::Tags::dictionary_reinit;


                // 3) send messages with local dofs to the right dict process
                for (const auto &rank_pair : buffers)
                  {
                    request.push_back(MPI_Request());
                    const int ierr = MPI_Isend(rank_pair.second.data(),
                                               rank_pair.second.size() * 2,
                                               DEAL_II_DOF_INDEX_MPI_TYPE,
                                               rank_pair.first,
                                               mpi_tag,
                                               comm,
                                               &request.back());
                    AssertThrowMPI(ierr);
                  }

                // 4) receive messages until all dofs in dict are processed
                while (this->local_size != dic_local_received)
                  {
                    // wait for an incoming message
                    MPI_Status status;
                    auto       ierr =
                      MPI_Probe(MPI_ANY_SOURCE, mpi_tag, comm, &status);
                    AssertThrowMPI(ierr);

                    // retrieve size of incoming message
                    int number_amount;
                    ierr = MPI_Get_count(&status,
                                         DEAL_II_DOF_INDEX_MPI_TYPE,
                                         &number_amount);
                    AssertThrowMPI(ierr);

                    const auto other_rank = status.MPI_SOURCE;
                    actually_owning_rank_list.push_back(other_rank);

                    // receive message
                    Assert(number_amount % 2 == 0, ExcInternalError());
                    std::vector<std::pair<types::global_dof_index,
                                          types::global_dof_index>>
                      buffer(number_amount / 2);
                    ierr = MPI_Recv(buffer.data(),
                                    number_amount,
                                    DEAL_II_DOF_INDEX_MPI_TYPE,
                                    status.MPI_SOURCE,
                                    status.MPI_TAG,
                                    comm,
                                    MPI_STATUS_IGNORE);
                    AssertThrowMPI(ierr);
                    // process message: loop over all intervals
                    for (auto interval : buffer)
                      {
#  ifdef DEBUG
                        for (types::global_dof_index i = interval.first;
                             i < interval.second;
                             i++)
                          Assert(actually_owning_ranks[i - local_range.first] ==
                                   numbers::invalid_unsigned_int,
                                 ExcInternalError());
                        Assert(interval.first >= local_range.first &&
                                 interval.first < local_range.second,
                               ExcInternalError());
                        Assert(interval.second > local_range.first &&
                                 interval.second <= local_range.second,
                               ExcInternalError());
#  endif

                        std::fill(actually_owning_ranks.data() +
                                    interval.first - local_range.first,
                                  actually_owning_ranks.data() +
                                    interval.second - local_range.first,
                                  other_rank);
                        dic_local_received += interval.second - interval.first;
                      }
                  }
              }
            else
              {
                // with gap: use a ConsensusAlgorithm to determine when all
                // dictionaries have been set up.

                // 3/4) use a ConsensusAlgorithm to send messages with local
                // dofs to the right dict process
                DictionaryPayLoad temp(buffers,
                                       actually_owning_ranks,
                                       local_range,
                                       actually_owning_rank_list);

                ConsensusAlgorithms::Selector<
                  std::pair<types::global_dof_index, types::global_dof_index>,
                  unsigned int>
                  consensus_algo(temp, comm);
                consensus_algo.run();
              }

            std::sort(actually_owning_rank_list.begin(),
                      actually_owning_rank_list.end());

            for (unsigned int i = 1; i < actually_owning_rank_list.size(); ++i)
              Assert(actually_owning_rank_list[i] >
                       actually_owning_rank_list[i - 1],
                     ExcInternalError());

            // 5) make sure that all messages have been sent
            if (request.size() > 0)
              {
                const int ierr = MPI_Waitall(request.size(),
                                             request.data(),
                                             MPI_STATUSES_IGNORE);
                AssertThrowMPI(ierr);
              }

#else
            (void)owned_indices;
            (void)comm;
#endif
          }

          /**
           * Translate a global dof index to the MPI rank in the dictionary
           * using `dofs_per_process`. We multiply by `stride_small_size` to
           * ensure a balance over the MPI ranks due to the grain size.
           */
          unsigned int
          dof_to_dict_rank(const types::global_dof_index i)
          {
            // note: this formula is also explicitly used in
            // get_index_offset(), so keep the two in sync
            return (i / dofs_per_process) * stride_small_size;
          }

          /**
           * Given an MPI rank id of an arbitrary processor, return the index
           * offset where the local range of that processor begins.
           */
          types::global_dof_index
          get_index_offset(const unsigned int rank)
          {
            return std::min(dofs_per_process *
                              static_cast<types::global_dof_index>(
                                (rank + stride_small_size - 1) /
                                stride_small_size),
                            size);
          }

          /**
           * Given the rank in the owned indices from `actually_owning_ranks`,
           * this returns the index of the rank in the
           * `actually_owning_rank_list`.
           */
          unsigned int
          get_owning_rank_index(const unsigned int rank_in_owned_indices,
                                const unsigned int guess = 0)
          {
            AssertIndexRange(guess, actually_owning_rank_list.size());
            if (actually_owning_rank_list[guess] == rank_in_owned_indices)
              return guess;
            else
              {
                auto it = std::lower_bound(actually_owning_rank_list.begin(),
                                           actually_owning_rank_list.end(),
                                           rank_in_owned_indices);
                Assert(it != actually_owning_rank_list.end(),
                       ExcInternalError());
                Assert(*it == rank_in_owned_indices, ExcInternalError());
                return it - actually_owning_rank_list.begin();
              }
          }

        private:
          /**
           * Compute the partition from the global size of the index space and
           * the number of ranks.
           */
          void
          partition(const IndexSet &owned_indices, const MPI_Comm &comm)
          {
#ifdef DEAL_II_WITH_MPI
            const unsigned int n_procs = n_mpi_processes(comm);
            const unsigned int my_rank = this_mpi_process(comm);

            size = owned_indices.size();

            Assert(size > 0, ExcNotImplemented());

            dofs_per_process = (size + n_procs - 1) / n_procs;
            if (dofs_per_process < range_minimum_grain_size)
              {
                dofs_per_process  = range_minimum_grain_size;
                stride_small_size = dofs_per_process * n_procs / size;
              }
            else
              stride_small_size = 1;
            local_range.first  = get_index_offset(my_rank);
            local_range.second = get_index_offset(my_rank + 1);

            local_size = local_range.second - local_range.first;

            actually_owning_ranks = {};
            actually_owning_ranks.resize(local_size,
                                         numbers::invalid_unsigned_int);
#else
            (void)owned_indices;
            (void)comm;
#endif
          }
        };



        /**
         * Specialization of ConsensusAlgorithms::Process for the context of
         * Utilities::MPI::compute_index_owner() and
         * Utilities::MPI::Partitioner::set_ghost_indices() with additional
         * payload.
         */
        class ConsensusAlgorithmsPayload
          : public ConsensusAlgorithms::Process<
              std::pair<types::global_dof_index, types::global_dof_index>,
              unsigned int>
        {
        public:
          /**
           * Constructor.
           */
          ConsensusAlgorithmsPayload(const IndexSet &owned_indices,
                                     const IndexSet &indices_to_look_up,
                                     const MPI_Comm &comm,
                                     std::vector<unsigned int> &owning_ranks,
                                     const bool track_index_requests = false)
            : owned_indices(owned_indices)
            , indices_to_look_up(indices_to_look_up)
            , comm(comm)
            , my_rank(this_mpi_process(comm))
            , n_procs(n_mpi_processes(comm))
            , track_index_requests(track_index_requests)
            , owning_ranks(owning_ranks)
          {
            dict.reinit(owned_indices, comm);
            requesters.resize(dict.actually_owning_rank_list.size());
          }

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
            std::vector<unsigned int> &request_buffer) override
          {
            unsigned int owner_index = 0;
            for (const auto &interval : buffer_recv)
              for (auto i = interval.first; i < interval.second; ++i)
                {
                  const unsigned int actual_owner =
                    dict.actually_owning_ranks[i - dict.local_range.first];
                  request_buffer.push_back(actual_owner);

                  if (track_index_requests)
                    append_index_origin(i, owner_index, other_rank);
                }
          }

          /**
           * Implementation of
           * Utilities::MPI::ConsensusAlgorithms::Process::compute_targets().
           */
          virtual std::vector<unsigned int>
          compute_targets() override
          {
            std::vector<unsigned int> targets;

            // 1) collect relevant processes and process local dict entries
            {
              unsigned int index       = 0;
              unsigned int owner_index = 0;
              for (auto i : indices_to_look_up)
                {
                  unsigned int other_rank = dict.dof_to_dict_rank(i);
                  if (other_rank == my_rank)
                    {
                      owning_ranks[index] =
                        dict.actually_owning_ranks[i - dict.local_range.first];
                      if (track_index_requests)
                        append_index_origin(i, owner_index, my_rank);
                    }
                  else if (targets.empty() || targets.back() != other_rank)
                    targets.push_back(other_rank);
                  index++;
                }
            }


            for (auto i : targets)
              {
                recv_indices[i]                    = {};
                indices_to_look_up_by_dict_rank[i] = {};
              }

            // 3) collect indices for each process
            {
              unsigned int index = 0;
              for (auto i : indices_to_look_up)
                {
                  unsigned int other_rank = dict.dof_to_dict_rank(i);
                  if (other_rank != my_rank)
                    {
                      recv_indices[other_rank].push_back(index);
                      indices_to_look_up_by_dict_rank[other_rank].push_back(i);
                    }
                  index++;
                }
            }

            Assert(targets.size() == recv_indices.size() &&
                     targets.size() == indices_to_look_up_by_dict_rank.size(),
                   ExcMessage("Size does not match!"));

            return targets;
          }

          /**
           * Implementation of
           * Utilities::MPI::ConsensusAlgorithms::Process::create_request().
           */
          virtual void
          create_request(const unsigned int other_rank,
                         std::vector<std::pair<types::global_dof_index,
                                               types::global_dof_index>>
                           &send_buffer) override
          {
            // create index set and compress data to be sent
            auto &   indices_i = indices_to_look_up_by_dict_rank[other_rank];
            IndexSet is(dict.size);
            is.add_indices(indices_i.begin(), indices_i.end());
            is.compress();

            for (auto interval = is.begin_intervals();
                 interval != is.end_intervals();
                 ++interval)
              send_buffer.emplace_back(*interval->begin(),
                                       interval->last() + 1);
          }

          /**
           * Implementation of
           * Utilities::MPI::ConsensusAlgorithms::Process::prepare_buffer_for_answer().
           */
          virtual void
          prepare_buffer_for_answer(
            const unsigned int         other_rank,
            std::vector<unsigned int> &recv_buffer) override
          {
            recv_buffer.resize(recv_indices[other_rank].size());
          }

          /**
           * Implementation of
           * Utilities::MPI::ConsensusAlgorithms::Process::read_answer().
           */
          virtual void
          read_answer(const unsigned int               other_rank,
                      const std::vector<unsigned int> &recv_buffer) override
          {
            Assert(recv_indices[other_rank].size() == recv_buffer.size(),
                   ExcMessage("Sizes do not match!"));

            for (unsigned int j = 0; j < recv_indices[other_rank].size(); j++)
              owning_ranks[recv_indices[other_rank][j]] = recv_buffer[j];
          }

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
          get_requesters()
          {
            Assert(track_index_requests,
                   ExcMessage("Must enable index range tracking in "
                              "constructor of ConsensusAlgorithmProcess"));

            std::map<unsigned int, dealii::IndexSet> requested_indices;

#ifdef DEAL_II_WITH_MPI

            static CollectiveMutex      mutex;
            CollectiveMutex::ScopedLock lock(mutex, comm);

            const int mpi_tag = Utilities::MPI::internal::Tags::
              consensus_algorithm_payload_get_requesters;

            // reserve enough slots for the requests ahead; depending on
            // whether the owning rank is one of the requesters or not, we
            // might have one less requests to execute, so fill the requests
            // on demand.
            std::vector<MPI_Request> send_requests;
            send_requests.reserve(requesters.size());

            // We use an integer vector for the data exchange. Since we send
            // data associated to intervals with different requesters, we will
            // need to send (a) the MPI rank of the requester, (b) the number
            // of intervals directed to this requester, and (c) a list of
            // intervals, i.e., two integers per interval. The number of items
            // sent in total can be deduced both via the MPI status message at
            // the receiver site as well as be counting the buckets from
            // different requesters.
            std::vector<std::vector<unsigned int>> send_data(requesters.size());
            for (unsigned int i = 0; i < requesters.size(); ++i)
              {
                // special code for our own indices
                if (dict.actually_owning_rank_list[i] == my_rank)
                  {
                    for (const auto &j : requesters[i])
                      {
                        const types::global_dof_index index_offset =
                          dict.get_index_offset(my_rank);
                        IndexSet &my_index_set = requested_indices[j.first];
                        my_index_set.set_size(owned_indices.size());
                        for (const auto &interval : j.second)
                          my_index_set.add_range(index_offset + interval.first,
                                                 index_offset +
                                                   interval.second);
                      }
                  }
                else
                  {
                    for (const auto &j : requesters[i])
                      {
                        send_data[i].push_back(j.first);
                        send_data[i].push_back(j.second.size());
                        for (const auto &interval : j.second)
                          {
                            send_data[i].push_back(interval.first);
                            send_data[i].push_back(interval.second);
                          }
                      }
                    send_requests.push_back(MPI_Request());
                    const int ierr =
                      MPI_Isend(send_data[i].data(),
                                send_data[i].size(),
                                MPI_UNSIGNED,
                                dict.actually_owning_rank_list[i],
                                mpi_tag,
                                comm,
                                &send_requests.back());
                    AssertThrowMPI(ierr);
                  }
              }

            // receive the data
            for (unsigned int c = 0; c < dict.n_dict_procs_in_owned_indices;
                 ++c)
              {
                // wait for an incoming message
                MPI_Status   status;
                unsigned int ierr =
                  MPI_Probe(MPI_ANY_SOURCE, mpi_tag, comm, &status);
                AssertThrowMPI(ierr);

                // retrieve size of incoming message
                int number_amount;
                ierr = MPI_Get_count(&status, MPI_UNSIGNED, &number_amount);
                AssertThrowMPI(ierr);

                // receive message
                Assert(number_amount % 2 == 0, ExcInternalError());
                std::vector<std::pair<unsigned int, unsigned int>> buffer(
                  number_amount / 2);
                ierr = MPI_Recv(buffer.data(),
                                number_amount,
                                MPI_UNSIGNED,
                                status.MPI_SOURCE,
                                status.MPI_TAG,
                                comm,
                                &status);
                AssertThrowMPI(ierr);

                // unpack the message and translate the dictionary-local
                // indices coming via MPI to the global index range
                const types::global_dof_index index_offset =
                  dict.get_index_offset(status.MPI_SOURCE);
                unsigned int offset = 0;
                while (offset < buffer.size())
                  {
                    AssertIndexRange(offset + buffer[offset].second,
                                     buffer.size());

                    IndexSet my_index_set(owned_indices.size());
                    for (unsigned int i = offset + 1;
                         i < offset + buffer[offset].second + 1;
                         ++i)
                      my_index_set.add_range(index_offset + buffer[i].first,
                                             index_offset + buffer[i].second);

                    // the underlying index set is able to merge ranges coming
                    // from different ranks due to the partitioning in the
                    // dictionary
                    IndexSet &index_set =
                      requested_indices[buffer[offset].first];
                    if (index_set.size() == 0)
                      index_set.set_size(owned_indices.size());
                    index_set.add_indices(my_index_set);

                    offset += buffer[offset].second + 1;
                  }
                AssertDimension(offset, buffer.size());
              }

            if (send_requests.size() > 0)
              {
                const auto ierr = MPI_Waitall(send_requests.size(),
                                              send_requests.data(),
                                              MPI_STATUSES_IGNORE);
                AssertThrowMPI(ierr);
              }


#  ifdef DEBUG
            for (const auto &it : requested_indices)
              {
                IndexSet copy_set = it.second;
                copy_set.subtract_set(owned_indices);
                Assert(copy_set.n_elements() == 0,
                       ExcInternalError(
                         "The indices requested from the current "
                         "MPI rank should be locally owned here!"));
              }
#  endif

#endif // DEAL_II_WITH_MPI

            return requested_indices;
          }

        private:
          /**
           * Stores the index request in the `requesters` field. We first find
           * out the owner of the index that was requested (using the guess in
           * `owner_index`, as we typically might look up on the same rank
           * several times in a row, which avoids the binary search in
           * Dictionary::get_owning_rank_index(). Once we know the rank of the
           * owner, we the vector entry with the rank of the request. Here, we
           * utilize the fact that requests are processed rank-by-rank, so we
           * can simply look at the end of the vector if there is already some
           * data stored or not. Finally, we build ranges, again using that
           * the index list is sorted and we therefore only need to append at
           * the end.
           */
          void
          append_index_origin(const types::global_dof_index index,
                              unsigned int &                owner_index,
                              const unsigned int            rank_of_request)
          {
            // remember who requested which index. We want to use an
            // std::vector with simple addressing, via a good guess from the
            // preceding index, rather than std::map, because this is an inner
            // loop and it avoids the map lookup in every iteration
            const unsigned int rank_of_owner =
              dict.actually_owning_ranks[index - dict.local_range.first];
            owner_index =
              dict.get_owning_rank_index(rank_of_owner, owner_index);
            if (requesters[owner_index].empty() ||
                requesters[owner_index].back().first != rank_of_request)
              requesters[owner_index].emplace_back(
                rank_of_request,
                std::vector<std::pair<unsigned int, unsigned int>>());
            if (requesters[owner_index].back().second.empty() ||
                requesters[owner_index].back().second.back().second !=
                  index - dict.local_range.first)
              requesters[owner_index].back().second.emplace_back(
                index - dict.local_range.first,
                index - dict.local_range.first + 1);
            else
              ++requesters[owner_index].back().second.back().second;
          }
        };

      } // namespace ComputeIndexOwner
    }   // namespace internal
  }     // namespace MPI
} // namespace Utilities

DEAL_II_NAMESPACE_CLOSE

#endif
