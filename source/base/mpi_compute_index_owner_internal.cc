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

#include <deal.II/base/config.h>

#include <deal.II/base/mpi.h>
#include <deal.II/base/mpi_compute_index_owner_internal.h>

#include <boost/serialization/utility.hpp>

DEAL_II_NAMESPACE_OPEN

namespace Utilities
{
  namespace MPI
  {
    namespace internal
    {
      namespace ComputeIndexOwner
      {
        const FlexibleIndexStorage::index_type
          FlexibleIndexStorage::invalid_index_value;



        FlexibleIndexStorage::FlexibleIndexStorage(const bool use_vector)
          : use_vector(use_vector)
          , size(0)
        {}



        void
        FlexibleIndexStorage::reinit(const bool        use_vector,
                                     const bool        index_range_contiguous,
                                     const std::size_t size)
        {
          this->use_vector = use_vector;
          this->size       = size;

          data = {};
          data_map.clear();

          // in case we have contiguous indices, only fill the vector upon
          // first request in `fill`
          if (!index_range_contiguous)
            data.resize(size, invalid_index_value);
        }



        void
        FlexibleIndexStorage::fill(
          const std::size_t                       start,
          const std::size_t                       end,
          const FlexibleIndexStorage::index_type &value)
        {
          AssertIndexRange(start, size);
          AssertIndexRange(end, size + 1);

          if (use_vector)
            {
              if (data.empty() && end > start)
                {
                  // in debug mode, we want to track whether we set all
                  // indices, so we first fill an invalid index and only later
                  // the actual ones, whereas we simply assign the given rank
                  // to the complete vector the first time we pass around in
                  // this function in release mode to avoid touching data
                  // unnecessarily (and overwrite the smaller pieces), as the
                  // locally owned part comes first
#ifdef DEBUG
                  data.resize(size, invalid_index_value);
                  std::fill(data.begin() + start, data.begin() + end, value);
#else
                  data.resize(size, value);
#endif
                }
              else
                {
                  AssertDimension(data.size(), size);
                  std::fill(data.begin() + start, data.begin() + end, value);
                }
            }
          else
            {
              for (auto i = start; i < end; ++i)
                data_map[i] = value;
            }
        }



        FlexibleIndexStorage::index_type &
        FlexibleIndexStorage::operator[](const std::size_t index)
        {
          AssertIndexRange(index, size);

          if (use_vector)
            {
              AssertDimension(data.size(), size);
              return data[index];
            }
          else
            {
              if (data_map.find(index) == data_map.end())
                data_map[index] = invalid_index_value;

              return data_map[index];
            }
        }



        FlexibleIndexStorage::index_type
        FlexibleIndexStorage::operator[](const std::size_t index) const
        {
          AssertIndexRange(index, size);

          if (use_vector)
            {
              AssertDimension(data.size(), size);
              return data[index];
            }
          else
            {
              if (data_map.find(index) == data_map.end())
                return invalid_index_value;

              return data_map.at(index);
            }
        }



        bool
        FlexibleIndexStorage::entry_has_been_set(const std::size_t index) const
        {
          AssertIndexRange(index, size);

          if (use_vector)
            {
              if (data.empty())
                return false;

              AssertDimension(data.size(), size);
              return data[index] != invalid_index_value;
            }
          else
            return data_map.find(index) != data_map.end();
        }



        DictionaryPayLoad::DictionaryPayLoad(
          const std::map<unsigned int,
                         std::vector<std::pair<types::global_dof_index,
                                               types::global_dof_index>>>
            &                   buffers,
          FlexibleIndexStorage &actually_owning_ranks,
          const std::pair<types::global_dof_index, types::global_dof_index>
            &                        local_range,
          std::vector<unsigned int> &actually_owning_rank_list)
          : buffers(buffers)
          , actually_owning_ranks(actually_owning_ranks)
          , local_range(local_range)
          , actually_owning_rank_list(actually_owning_rank_list)
        {
          Assert(local_range.first <= local_range.second, ExcInternalError());
        }



        std::vector<unsigned int>
        DictionaryPayLoad::compute_targets()
        {
          std::vector<unsigned int> targets;
          for (const auto &rank_pair : buffers)
            targets.push_back(rank_pair.first);

          return targets;
        }



        void
        DictionaryPayLoad::create_request(
          const unsigned int                               other_rank,
          std::vector<std::pair<types::global_dof_index,
                                types::global_dof_index>> &send_buffer)
        {
          send_buffer = this->buffers.at(other_rank);
        }



        void
        DictionaryPayLoad::answer_request(
          const unsigned int                                     other_rank,
          const std::vector<std::pair<types::global_dof_index,
                                      types::global_dof_index>> &buffer_recv,
          std::vector<unsigned int> &                            request_buffer)
        {
          (void)request_buffer; // not needed


          // process message: loop over all intervals
          for (auto interval : buffer_recv)
            {
#ifdef DEBUG
              for (types::global_dof_index i = interval.first;
                   i < interval.second;
                   i++)
                Assert(
                  actually_owning_ranks.entry_has_been_set(
                    i - local_range.first) == false,
                  ExcMessage(
                    "Multiple processes seem to own the same global index. "
                    "A possible reason is that the sets of locally owned "
                    "indices are not distinct."));
              Assert(interval.first < interval.second, ExcInternalError());
              Assert(
                local_range.first <= interval.first &&
                  interval.second <= local_range.second,
                ExcMessage(
                  "The specified interval is not handled by the current process."));
#endif
              actually_owning_ranks.fill(interval.first - local_range.first,
                                         interval.second - local_range.first,
                                         other_rank);
            }
          actually_owning_rank_list.push_back(other_rank);
        }



        void
        Dictionary::reinit(const IndexSet &owned_indices, const MPI_Comm &comm)
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

          const auto owned_indices_size_actual =
            Utilities::MPI::sum(owned_indices.n_elements(), comm);

          actually_owning_ranks.reinit((owned_indices_size_actual *
                                        sparsity_factor) > owned_indices.size(),
                                       owned_indices_size_actual ==
                                         owned_indices.size(),
                                       locally_owned_size);

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

              AssertThrow(index_range.second <= size, ExcInternalError());

              while (index_range.first != index_range.second)
                {
                  Assert(index_range.first < index_range.second,
                         ExcInternalError());

                  const unsigned int owner =
                    dof_to_dict_rank(index_range.first);

                  // this explicitly picks up the formula of
                  // dof_to_dict_rank, so the two places must be in sync
                  const types::global_dof_index next_index =
                    std::min(get_index_offset(owner + 1), index_range.second);

                  Assert(next_index > index_range.first, ExcInternalError());

#  ifdef DEBUG
                  // make sure that the owner is the same on the current
                  // interval
                  for (types::global_dof_index i = index_range.first + 1;
                       i < next_index;
                       ++i)
                    AssertDimension(owner, dof_to_dict_rank(i));
#  endif

                  // add the interval, either to the local range or into a
                  // buffer to be sent to another processor
                  if (owner == my_rank)
                    {
                      actually_owning_ranks.fill(index_range.first -
                                                   local_range.first,
                                                 next_index - local_range.first,
                                                 my_rank);
                      dic_local_received += next_index - index_range.first;
                      if (actually_owning_rank_list.empty())
                        actually_owning_rank_list.push_back(my_rank);
                    }
                  else
                    buffers[owner].emplace_back(index_range.first, next_index);

                  index_range.first = next_index;
                }
            }

          n_dict_procs_in_owned_indices = buffers.size();
          std::vector<MPI_Request> request;

          // Check if index set space is partitioned globally without gaps.
          if (owned_indices_size_actual == owned_indices.size())
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
              while (this->locally_owned_size != dic_local_received)
                {
                  // wait for an incoming message
                  MPI_Status status;
                  int ierr = MPI_Probe(MPI_ANY_SOURCE, mpi_tag, comm, &status);
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
                        Assert(actually_owning_ranks.entry_has_been_set(
                                 i - local_range.first) == false,
                               ExcInternalError());
                      Assert(interval.first >= local_range.first &&
                               interval.first < local_range.second,
                             ExcInternalError());
                      Assert(interval.second > local_range.first &&
                               interval.second <= local_range.second,
                             ExcInternalError());
#  endif

                      actually_owning_ranks.fill(interval.first -
                                                   local_range.first,
                                                 interval.second -
                                                   local_range.first,
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
                std::vector<
                  std::pair<types::global_dof_index, types::global_dof_index>>,
                std::vector<unsigned int>>
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



        void
        Dictionary::partition(const IndexSet &owned_indices,
                              const MPI_Comm &comm)
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

          locally_owned_size = local_range.second - local_range.first;
#else
          (void)owned_indices;
          (void)comm;
#endif
        }


        ConsensusAlgorithmsPayload::ConsensusAlgorithmsPayload(
          const IndexSet &           owned_indices,
          const IndexSet &           indices_to_look_up,
          const MPI_Comm &           comm,
          std::vector<unsigned int> &owning_ranks,
          const bool                 track_index_requests)
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



        void
        ConsensusAlgorithmsPayload::answer_request(
          const unsigned int                                     other_rank,
          const std::vector<std::pair<types::global_dof_index,
                                      types::global_dof_index>> &buffer_recv,
          std::vector<unsigned int> &                            request_buffer)
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



        std::vector<unsigned int>
        ConsensusAlgorithmsPayload::compute_targets()
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



        void
        ConsensusAlgorithmsPayload::create_request(
          const unsigned int                               other_rank,
          std::vector<std::pair<types::global_dof_index,
                                types::global_dof_index>> &send_buffer)
        {
          // create index set and compress data to be sent
          auto &   indices_i = indices_to_look_up_by_dict_rank[other_rank];
          IndexSet is(dict.size);
          is.add_indices(indices_i.begin(), indices_i.end());
          is.compress();

          for (auto interval = is.begin_intervals();
               interval != is.end_intervals();
               ++interval)
            send_buffer.emplace_back(*interval->begin(), interval->last() + 1);
        }



        void
        ConsensusAlgorithmsPayload::read_answer(
          const unsigned int               other_rank,
          const std::vector<unsigned int> &recv_buffer)
        {
          Assert(recv_buffer.size() == recv_indices[other_rank].size(),
                 ExcInternalError());
          Assert(recv_indices[other_rank].size() == recv_buffer.size(),
                 ExcMessage("Sizes do not match!"));

          for (unsigned int j = 0; j < recv_indices[other_rank].size(); ++j)
            owning_ranks[recv_indices[other_rank][j]] = recv_buffer[j];
        }



        std::map<unsigned int, IndexSet>
        ConsensusAlgorithmsPayload::get_requesters()
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
                                               index_offset + interval.second);
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
                  const int ierr = MPI_Isend(send_data[i].data(),
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
          for (unsigned int c = 0; c < dict.n_dict_procs_in_owned_indices; ++c)
            {
              // wait for an incoming message
              MPI_Status status;
              int ierr = MPI_Probe(MPI_ANY_SOURCE, mpi_tag, comm, &status);
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
                  IndexSet &index_set = requested_indices[buffer[offset].first];
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



        void
        ConsensusAlgorithmsPayload::append_index_origin(
          const types::global_dof_index index,
          unsigned int &                owner_index,
          const unsigned int            rank_of_request)
        {
          // remember who requested which index. We want to use an
          // std::vector with simple addressing, via a good guess from the
          // preceding index, rather than std::map, because this is an inner
          // loop and it avoids the map lookup in every iteration
          const unsigned int rank_of_owner =
            dict.actually_owning_ranks[index - dict.local_range.first];
          owner_index = dict.get_owning_rank_index(rank_of_owner, owner_index);
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
      } // namespace ComputeIndexOwner
    }   // namespace internal
  }     // namespace MPI
} // namespace Utilities

DEAL_II_NAMESPACE_CLOSE
