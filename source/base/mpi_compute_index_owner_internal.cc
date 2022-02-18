// ---------------------------------------------------------------------
//
// Copyright (C) 2019 - 2021 by the deal.II authors
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
        DictionaryPayLoad::DictionaryPayLoad(
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
                  actually_owning_ranks[i - local_range.first] ==
                    numbers::invalid_unsigned_int,
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
              std::fill(actually_owning_ranks.data() + interval.first -
                          local_range.first,
                        actually_owning_ranks.data() + interval.second -
                          local_range.first,
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

                      std::fill(actually_owning_ranks.data() + interval.first -
                                  local_range.first,
                                actually_owning_ranks.data() + interval.second -
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

          actually_owning_ranks = {};
          actually_owning_ranks.resize(locally_owned_size,
                                       numbers::invalid_unsigned_int);
#else
          (void)owned_indices;
          (void)comm;
#endif
        }


      } // namespace ComputeIndexOwner
    }   // namespace internal
  }     // namespace MPI
} // namespace Utilities

DEAL_II_NAMESPACE_CLOSE
