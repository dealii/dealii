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
       * An internal namespace used for Utilities::MPI::compute_index_owner().
       */
      namespace ComputeIndexOwner
      {
        struct Dictionary
        {
          static const unsigned int tag_setup = 11;

          std::vector<unsigned int> actually_owning_ranks;

          types::global_dof_index dofs_per_process;
          std::pair<types::global_dof_index, types::global_dof_index>
                                  local_range;
          types::global_dof_index local_size;
          types::global_dof_index size;

          void
          reinit(const IndexSet &owned_indices, const MPI_Comm &comm)
          {
            this->partition(owned_indices, comm);

#ifdef DEAL_II_WITH_MPI
            unsigned int my_rank = this_mpi_process(comm);

            types::global_dof_index              dic_local_rececived = 0;
            std::map<unsigned int, unsigned int> relevant_procs_map;

            // 2) collect relevant processes and process local dict entries
            {
              std::vector<unsigned int> relevant_procs;
              for (auto i : owned_indices)
                {
                  unsigned int other_rank = this->dof_to_dict_rank(i);
                  if (other_rank == my_rank)
                    {
                      this->actually_owning_ranks[i - this->local_range.first] =
                        my_rank;
                      dic_local_rececived++;
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

            const unsigned int n_relevant_procs = relevant_procs_map.size();
            std::vector<std::vector<
              std::pair<types::global_dof_index, types::global_dof_index>>>
                                     buffers(n_relevant_procs);
            std::vector<MPI_Request> request(n_relevant_procs);

            // 3) send messages with local dofs to the right dict process
            {
              std::vector<std::vector<types::global_dof_index>> temp(
                n_relevant_procs);

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
            while (this->local_size != dic_local_rececived)
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
                      dic_local_rececived++;
                    }
              }

            // 5) make sure that all messages have been sent
            const auto ierr = MPI_Waitall(n_relevant_procs,
                                          request.data(),
                                          MPI_STATUSES_IGNORE);
            AssertThrowMPI(ierr);
#else
            (void)owned_indices;
            (void)comm;
#endif
          }

          unsigned int
          dof_to_dict_rank(const types::global_dof_index i)
          {
            return i / dofs_per_process;
          }

        private:
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

        class ConsensusAlgorithmProcess
          : public dealii::Utilities::MPI::
              ConsensusAlgorithmProcess<types::global_dof_index, unsigned int>
        {
        public:
          ConsensusAlgorithmProcess(const IndexSet &owned_indices,
                                    const IndexSet &indices_to_look_up,
                                    const MPI_Comm &comm,
                                    std::vector<unsigned int> &owning_ranks)
            : owned_indices(owned_indices)
            , indices_to_look_up(indices_to_look_up)
            , comm(comm)
            , my_rank(this_mpi_process(comm))
            , n_procs(n_mpi_processes(comm))
            , owning_ranks(owning_ranks)
          {
            this->dict.reinit(owned_indices, comm);
          }

          const IndexSet &           owned_indices;
          const IndexSet &           indices_to_look_up;
          const MPI_Comm &           comm;
          const unsigned int         my_rank;
          const unsigned int         n_procs;
          std::vector<unsigned int> &owning_ranks;

          Dictionary dict;

          std::map<unsigned int, std::vector<types::global_dof_index>> temp;
          std::map<unsigned int, std::vector<unsigned int>> recv_indices;

          virtual void
          process_request(
            const unsigned int                          other_rank,
            const std::vector<types::global_dof_index> &buffer_recv,
            std::vector<unsigned int> &                 request_buffer) override
          {
            (void)other_rank;
            Assert(buffer_recv.size() % 2 == 0, ExcInternalError());
            for (unsigned int j = 0; j < buffer_recv.size(); j += 2)
              for (auto i = buffer_recv[j]; i < buffer_recv[j + 1]; i++)
                request_buffer.push_back(
                  dict.actually_owning_ranks[i - dict.local_range.first]);
          }

          virtual std::vector<unsigned int>
          compute_targets() override
          {
            std::vector<unsigned int> targets;

            // 1) collect relevant processes and process local dict entries
            {
              unsigned int index = 0;
              for (auto i : indices_to_look_up)
                {
                  unsigned int other_rank = dict.dof_to_dict_rank(i);
                  if (other_rank == my_rank)
                    owning_ranks[index] =
                      dict.actually_owning_ranks[i - dict.local_range.first];
                  else if (targets.empty() || targets.back() != other_rank)
                    targets.push_back(other_rank);
                  index++;
                }
            }


            for (auto i : targets)
              {
                recv_indices[i] = {};
                temp[i]         = {};
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
                      temp[other_rank].push_back(i);
                    }
                  index++;
                }
            }

            Assert(targets.size() == recv_indices.size() &&
                     targets.size() == temp.size(),
                   ExcMessage("Size does not match!"));

            return targets;
          }

          virtual void
          pack_recv_buffer(
            const int                             other_rank,
            std::vector<types::global_dof_index> &send_buffer) override
          {
            // create index set and compress data to be sent
            auto &   indices_i = temp[other_rank];
            IndexSet is(dict.size);
            is.add_indices(indices_i.begin(), indices_i.end());
            is.compress();

            for (auto interval = is.begin_intervals();
                 interval != is.end_intervals();
                 interval++)
              {
                send_buffer.push_back(*interval->begin());
                send_buffer.push_back(interval->last() + 1);
              }
          }

          virtual void
          prepare_recv_buffer(const int                  other_rank,
                              std::vector<unsigned int> &recv_buffer) override
          {
            recv_buffer.resize(recv_indices[other_rank].size());
          }

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
        };

      } // namespace ComputeIndexOwner
    }   // namespace internal
  }     // namespace MPI
} // namespace Utilities

DEAL_II_NAMESPACE_CLOSE

#endif
