// ---------------------------------------------------------------------
//
// Copyright (C) 2019 by the deal.II authors
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
         * Dictionary class with basic partitioning in terms of a single
         * interval of fixed size known to all MPI ranks for two-stage index
         * lookup.
         */
        struct Dictionary
        {
          /**
           * A tag attached to the MPI communication during the dictionary
           * lookup
           */
          static const unsigned int tag_setup = 11;

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
           * the dictionary, computed from `dofs_per_process` and the current
           * MPI rank.
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
           * Set up the dictionary by computing the partitioning from the
           * global size and sending the rank information on locally owned
           * ranges to the owner of the dictionary part.
           */
          void
          reinit(const IndexSet &owned_indices, const MPI_Comm &comm)
          {
            this->partition(owned_indices, comm);

#ifdef DEAL_II_WITH_MPI
            unsigned int my_rank = this_mpi_process(comm);

            types::global_dof_index              dic_local_received = 0;
            std::map<unsigned int, unsigned int> relevant_procs_map;

            // 2) collect relevant processes and process local dict entries
            {
              std::vector<unsigned int> relevant_procs;
              for (auto i : owned_indices)
                {
                  unsigned int other_rank = this->dof_to_dict_rank(i);
                  if (other_rank == my_rank)
                    {
                      actually_owning_ranks[i - local_range.first] = my_rank;
                      dic_local_received++;
                      if (actually_owning_rank_list.empty())
                        actually_owning_rank_list.push_back(my_rank);
                    }
                  else if (relevant_procs.empty() ||
                           relevant_procs.back() != other_rank)
                    relevant_procs.push_back(other_rank);
                }

              {
                unsigned int c = 0;
                for (auto i : relevant_procs)
                  relevant_procs_map[i] = c++;
              }
            }

            n_dict_procs_in_owned_indices = relevant_procs_map.size();
            std::vector<std::vector<
              std::pair<types::global_dof_index, types::global_dof_index>>>
                                     buffers(n_dict_procs_in_owned_indices);
            std::vector<MPI_Request> request(n_dict_procs_in_owned_indices);

            // 3) send messages with local dofs to the right dict process
            {
              std::vector<std::vector<types::global_dof_index>> temp(
                n_dict_procs_in_owned_indices);

              // collect dofs of each dict process
              for (auto i : owned_indices)
                {
                  unsigned int other_rank = this->dof_to_dict_rank(i);
                  if (other_rank != my_rank)
                    temp[relevant_procs_map[other_rank]].push_back(i);
                }

              // send dofs to each process
              for (auto rank_pair : relevant_procs_map)
                {
                  const int rank  = rank_pair.first;
                  const int index = rank_pair.second;

                  // create index set and compress data to be sent
                  auto &   indices_i = temp[index];
                  IndexSet is(this->size);
                  is.add_indices(indices_i.begin(), indices_i.end());
                  is.compress();

                  // translate index set to a list of pairs
                  auto &buffer = buffers[index];
                  for (auto interval = is.begin_intervals();
                       interval != is.end_intervals();
                       interval++)
                    buffer.emplace_back(*interval->begin(),
                                        interval->last() + 1);

                  // send data
                  const auto ierr = MPI_Isend(buffer.data(),
                                              buffer.size() * 2,
                                              DEAL_II_DOF_INDEX_MPI_TYPE,
                                              rank,
                                              tag_setup,
                                              comm,
                                              &request[index]);
                  AssertThrowMPI(ierr);
                }
            }

            // 4) receive messages until all dofs in dict are processed
            while (this->local_size != dic_local_received)
              {
                // wait for an incoming message
                MPI_Status status;
                auto ierr = MPI_Probe(MPI_ANY_SOURCE, tag_setup, comm, &status);
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
                std::vector<
                  std::pair<types::global_dof_index, types::global_dof_index>>
                  buffer(number_amount / 2);
                ierr = MPI_Recv(buffer.data(),
                                number_amount,
                                DEAL_II_DOF_INDEX_MPI_TYPE,
                                other_rank,
                                tag_setup,
                                comm,
                                &status);
                AssertThrowMPI(ierr);

                // process message: loop over all intervals
                for (auto interval : buffer)
                  for (types::global_dof_index i = interval.first;
                       i < interval.second;
                       i++)
                    {
                      this->actually_owning_ranks[i - this->local_range.first] =
                        other_rank;
                      dic_local_received++;
                    }
              }

            std::sort(actually_owning_rank_list.begin(),
                      actually_owning_rank_list.end());

            for (unsigned int i = 1; i < actually_owning_rank_list.size(); ++i)
              Assert(actually_owning_rank_list[i] >
                       actually_owning_rank_list[i - 1],
                     ExcInternalError());

            // 5) make sure that all messages have been sent
            const auto ierr =
              MPI_Waitall(request.size(), request.data(), MPI_STATUSES_IGNORE);
            AssertThrowMPI(ierr);
#else
            (void)owned_indices;
            (void)comm;
#endif
          }

          /**
           * Translate a global dof index to the MPI rank in the dictionary
           * using `dofs_per_process`.
           */
          unsigned int
          dof_to_dict_rank(const types::global_dof_index i)
          {
            return i / dofs_per_process;
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

            size              = owned_indices.size();
            dofs_per_process  = (size + n_procs - 1) / n_procs;
            local_range.first = std::min(dofs_per_process * my_rank, size);
            local_range.second =
              std::min(dofs_per_process * (my_rank + 1), size);
            local_size = local_range.second - local_range.first;

            actually_owning_ranks.resize(local_size);
#else
            (void)owned_indices;
            (void)comm;
#endif
          }
        };



        /**
         * Specialization of ConsensusAlgorithmProcess for the context of
         * Utilities::MPI::compute_index_owner() and
         * Utilities::MPI::Partitioner::set_ghost_indices() with additional
         * payload.
         */
        class ConsensusAlgorithmPayload
          : public ConsensusAlgorithmProcess<
              std::pair<types::global_dof_index, types::global_dof_index>,
              unsigned int>
        {
        public:
          /**
           * Constructor.
           */
          ConsensusAlgorithmPayload(const IndexSet &owned_indices,
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
           * Utilities::MPI::ConsensusAlgorithmProcess::process_request(),
           * adding the owner of a particular index in request_buffer (and
           * keeping track of who requested a particular index in case that
           * information is also desired).
           */
          virtual void
          process_request(
            const unsigned int                                     other_rank,
            const std::vector<std::pair<types::global_dof_index,
                                        types::global_dof_index>> &buffer_recv,
            std::vector<unsigned int> &request_buffer) override
          {
            unsigned int owner_index = 0;
            for (const auto interval : buffer_recv)
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
           * Utilities::MPI::ConsensusAlgorithmProcess::compute_targets().
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
           * Utilities::MPI::ConsensusAlgorithmProcess::pack_recv_buffer().
           */
          virtual void
          pack_recv_buffer(const int other_rank,
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
                 interval++)
              send_buffer.emplace_back(*interval->begin(),
                                       interval->last() + 1);
          }

          /**
           * Implementation of
           * Utilities::MPI::ConsensusAlgorithmProcess::prepare_recv_buffer().
           */
          virtual void
          prepare_recv_buffer(const int                  other_rank,
                              std::vector<unsigned int> &recv_buffer) override
          {
            recv_buffer.resize(recv_indices[other_rank].size());
          }

          /**
           * Implementation of
           * Utilities::MPI::ConsensusAlgorithmProcess::unpack_recv_buffer().
           */
          virtual void
          unpack_recv_buffer(
            const int                        other_rank,
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
                   ExcMessage("Must enable index range tracking in"
                              "constructor of ConsensusAlgorithmProcess"));

            std::map<unsigned int, dealii::IndexSet> requested_indices;

#ifdef DEAL_II_WITH_MPI

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
                          static_cast<types::global_dof_index>(my_rank) *
                          dict.dofs_per_process;
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
                                1021,
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
                  MPI_Probe(MPI_ANY_SOURCE, 1021, comm, &status);
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
                                1021,
                                comm,
                                &status);
                AssertThrowMPI(ierr);

                // unpack the message and translate the dictionary-local
                // indices coming via MPI to the global index range
                const types::global_dof_index index_offset =
                  static_cast<types::global_dof_index>(status.MPI_SOURCE) *
                  dict.dofs_per_process;
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
              MPI_Waitall(send_requests.size(),
                          send_requests.data(),
                          MPI_STATUSES_IGNORE);

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
