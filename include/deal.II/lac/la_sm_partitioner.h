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

#ifndef dealii_la_sm_partitioner_h
#define dealii_la_sm_partitioner_h


#include <deal.II/base/config.h>

#include <deal.II/base/mpi_compute_index_owner_internal.h>

#include <deal.II/lac/communication_pattern_base.h>

DEAL_II_NAMESPACE_OPEN


namespace LinearAlgebra
{
  namespace SharedMPI
  {
    template <typename Number = double>
    class Partitioner : public LinearAlgebra::CommunicationPatternBase
    {
    public:
      Partitioner(const MPI_Comm &comm,
                  const MPI_Comm &comm_sm,
                  const IndexSet &is_locally_owned,
                  const IndexSet &is_locally_ghost)
        : comm(comm)
        , comm_sm(comm_sm)
      {
        reinit(is_locally_owned, is_locally_ghost);
      }

      const MPI_Comm &
      get_mpi_communicator() const override
      {
        return comm;
      }

      const MPI_Comm &
      get_sm_mpi_communicator() const
      {
        return comm_sm;
      }

      void
      reinit(const IndexSet &is_locally_owned,
             const IndexSet &is_locally_ghost,
             const MPI_Comm &communicator) override
      {
        Assert(false, ExcNotImplemented());
        (void)is_locally_owned;
        (void)is_locally_ghost;
        (void)communicator;
      }

      void
      reinit(const IndexSet &is_locally_owned, const IndexSet &is_locally_ghost)
      {
        this->n_local_elements = is_locally_owned.n_elements();

        std::vector<unsigned int> sm_ranks(
          Utilities::MPI::n_mpi_processes(comm_sm));

        const unsigned int rank = Utilities::MPI::this_mpi_process(comm);

        MPI_Allgather(
          &rank, 1, MPI_UNSIGNED, sm_ranks.data(), 1, MPI_UNSIGNED, comm_sm);

        std::vector<unsigned int> owning_ranks_of_ghosts(
          is_locally_ghost.n_elements());

        Utilities::MPI::internal::ComputeIndexOwner::ConsensusAlgorithmsPayload
          process(is_locally_owned,
                  is_locally_ghost,
                  comm,
                  owning_ranks_of_ghosts,
                  true);

        Utilities::MPI::ConsensusAlgorithms::Selector<
          std::pair<types::global_dof_index, types::global_dof_index>,
          unsigned int>
          consensus_algorithm(process, comm);
        consensus_algorithm.run();

        {
          std::map<unsigned int, std::vector<types::global_dof_index>>
            rank_to_local_indices;

          for (unsigned int i = 0; i < owning_ranks_of_ghosts.size(); i++)
            rank_to_local_indices[owning_ranks_of_ghosts[i]].push_back(i);

          unsigned int offset = 0;


          for (const auto &rank_and_local_indices : rank_to_local_indices)
            {
              const auto ptr = std::find(sm_ranks.begin(),
                                         sm_ranks.end(),
                                         rank_and_local_indices.first);

              if (ptr == sm_ranks.end())
                {
                  // remote process
                  recv_remote_ranks.push_back(rank_and_local_indices.first);
                  recv_remote_ptr.push_back(
                    recv_remote_ptr.back() +
                    rank_and_local_indices.second.size());
                }
              else
                {
                  // shared process
                  recv_sm_ranks.push_back(std::distance(sm_ranks.begin(), ptr));
                  recv_sm_ptr.push_back(recv_sm_ptr.back() +
                                        rank_and_local_indices.second.size());
                  recv_sm_offset.push_back(is_locally_owned.n_elements() +
                                           offset);
                }
              offset += rank_and_local_indices.second.size();
            }
          recv_remote_req.resize(recv_remote_ranks.size());
          recv_sm_req.resize(recv_sm_ranks.size());

          recv_sm_indices.resize(recv_sm_ptr.back());
        }

        {
          const auto rank_to_global_indices = process.get_requesters();

          for (const auto &rank_and_global_indices : rank_to_global_indices)
            {
              const auto ptr = std::find(sm_ranks.begin(),
                                         sm_ranks.end(),
                                         rank_and_global_indices.first);

              if (ptr == sm_ranks.end())
                {
                  // remote process
                  send_remote_ranks.push_back(rank_and_global_indices.first);

                  for (const auto &i : rank_and_global_indices.second)
                    send_remote_indices.push_back(
                      is_locally_owned.index_within_set(i));

                  send_remote_ptr.push_back(send_remote_indices.size());
                }
              else
                {
                  // shared process
                  send_sm_ranks.push_back(std::distance(sm_ranks.begin(), ptr));

                  for (const auto &i : rank_and_global_indices.second)
                    send_sm_indices.push_back(
                      is_locally_owned.index_within_set(i));

                  send_sm_ptr.push_back(send_sm_indices.size());
                }
            }
          send_remote_req.resize(send_remote_ranks.size());
          send_sm_req.resize(send_sm_ranks.size());
        }

        {
          for (unsigned int i = 0; i < send_sm_ranks.size(); i++)
            MPI_Isend(send_sm_indices.data() + send_sm_ptr[i],
                      send_sm_ptr[i + 1] - send_sm_ptr[i],
                      MPI_UNSIGNED,
                      send_sm_ranks[i],
                      2,
                      comm_sm,
                      send_sm_req.data() + i);

          for (unsigned int i = 0; i < recv_sm_ranks.size(); i++)
            MPI_Irecv(recv_sm_indices.data() + recv_sm_ptr[i],
                      recv_sm_ptr[i + 1] - recv_sm_ptr[i],
                      MPI_UNSIGNED,
                      recv_sm_ranks[i],
                      2,
                      comm_sm,
                      recv_sm_req.data() + i);

          MPI_Waitall(recv_sm_req.size(),
                      recv_sm_req.data(),
                      MPI_STATUSES_IGNORE);
          MPI_Waitall(send_sm_req.size(),
                      send_sm_req.data(),
                      MPI_STATUSES_IGNORE);
        }

        {
          send_sm_offset.resize(send_sm_ranks.size());

          for (unsigned int i = 0; i < send_sm_ranks.size(); i++)
            MPI_Irecv(send_sm_offset.data() + i,
                      1,
                      MPI_UNSIGNED,
                      send_sm_ranks[i],
                      3,
                      comm_sm,
                      send_sm_req.data() + i);

          for (unsigned int i = 0; i < recv_sm_ranks.size(); i++)
            MPI_Isend(recv_sm_offset.data() + i,
                      1,
                      MPI_UNSIGNED,
                      recv_sm_ranks[i],
                      3,
                      comm_sm,
                      recv_sm_req.data() + i);

          MPI_Waitall(recv_sm_req.size(),
                      recv_sm_req.data(),
                      MPI_STATUSES_IGNORE);
          MPI_Waitall(send_sm_req.size(),
                      send_sm_req.data(),
                      MPI_STATUSES_IGNORE);
        }

        buffer.resize(send_remote_ptr.back());
      }

      void
      update_ghost_values_start(
        Number *               data_this,
        std::vector<Number *> &data_others,
        const unsigned int     communication_channel = 0) const
      {
        (void)data_others;

        int dummy;
        for (unsigned int i = 0; i < recv_sm_ranks.size(); i++)
          MPI_Irecv(&dummy,
                    0,
                    MPI_INT,
                    recv_sm_ranks[i],
                    communication_channel + 1,
                    comm_sm,
                    recv_sm_req.data() + i);

        for (unsigned int i = 0; i < send_sm_ranks.size(); i++)
          MPI_Isend(&dummy,
                    0,
                    MPI_INT,
                    send_sm_ranks[i],
                    communication_channel + 1,
                    comm_sm,
                    send_sm_req.data() + i);

        for (unsigned int i = 0; i < recv_remote_ranks.size(); i++)
          MPI_Irecv(data_this + recv_remote_ptr[i] + n_local_elements,
                    recv_remote_ptr[i + 1] - recv_remote_ptr[i],
                    Utilities::MPI::internal::mpi_type_id(data_this),
                    recv_remote_ranks[i],
                    communication_channel + 0,
                    comm,
                    recv_remote_req.data() + i);

        for (unsigned int i = 0; i < send_remote_ranks.size(); i++)
          {
            for (unsigned int j = send_remote_ptr[i];
                 j < send_remote_ptr[i + 1];
                 j++)
              buffer[j] = data_this[send_remote_indices[j]];

            MPI_Isend(buffer.data() + send_remote_ptr[i],
                      send_remote_ptr[i + 1] - send_remote_ptr[i],
                      Utilities::MPI::internal::mpi_type_id(data_this),
                      send_remote_ranks[i],
                      communication_channel + 0,
                      comm,
                      send_remote_req.data() + i);
          }
      }

      void
      update_ghost_values_finish(Number *               data_this,
                                 std::vector<Number *> &data_others) const
      {
        for (unsigned int c = 0; c < recv_sm_ranks.size(); c++)
          {
            int i;
            MPI_Waitany(recv_sm_req.size(),
                        recv_sm_req.data(),
                        &i,
                        MPI_STATUS_IGNORE);

            const Number *__restrict__ data_others_ptr =
              data_others[recv_sm_ranks[i]];
            Number *__restrict__ data_this_ptr = data_this;

            for (unsigned int j = recv_sm_ptr[i], k = recv_sm_offset[i];
                 j < recv_sm_ptr[i + 1];
                 j++, k++)
              data_this_ptr[k] = data_others_ptr[recv_sm_indices[j]];
          }

        MPI_Waitall(send_sm_req.size(),
                    send_sm_req.data(),
                    MPI_STATUSES_IGNORE);
        MPI_Waitall(send_remote_req.size(),
                    send_remote_req.data(),
                    MPI_STATUSES_IGNORE);
        MPI_Waitall(recv_remote_req.size(),
                    recv_remote_req.data(),
                    MPI_STATUSES_IGNORE);
      }

      void
      update_ghost_values(Number *               data_this,
                          std::vector<Number *> &data_others) const
      {
        update_ghost_values_start(data_this, data_others);
        update_ghost_values_finish(data_this, data_others);
      }

      void
      compress_start(Number *               data_this,
                     std::vector<Number *> &data_others,
                     const unsigned int     communication_channel = 0) const
      {
        (void)data_others;

        int dummy;
        for (unsigned int i = 0; i < recv_sm_ranks.size(); i++)
          MPI_Isend(&dummy,
                    0,
                    MPI_INT,
                    recv_sm_ranks[i],
                    communication_channel + 1,
                    comm_sm,
                    recv_sm_req.data() + i);

        for (unsigned int i = 0; i < send_sm_ranks.size(); i++)
          MPI_Irecv(&dummy,
                    0,
                    MPI_INT,
                    send_sm_ranks[i],
                    communication_channel + 1,
                    comm_sm,
                    send_sm_req.data() + i);

        for (unsigned int i = 0; i < recv_remote_ranks.size(); i++)
          MPI_Isend(data_this + recv_remote_ptr[i] + n_local_elements,
                    recv_remote_ptr[i + 1] - recv_remote_ptr[i],
                    Utilities::MPI::internal::mpi_type_id(data_this),
                    recv_remote_ranks[i],
                    communication_channel + 0,
                    comm,
                    recv_remote_req.data() + i);

        for (unsigned int i = 0; i < send_remote_ranks.size(); i++)
          MPI_Irecv(buffer.data() + send_remote_ptr[i],
                    send_remote_ptr[i + 1] - send_remote_ptr[i],
                    Utilities::MPI::internal::mpi_type_id(data_this),
                    send_remote_ranks[i],
                    communication_channel + 0,
                    comm,
                    send_remote_req.data() + i);
      }

      void
      compress_finish(Number *               data_this,
                      std::vector<Number *> &data_others) const
      {
        for (unsigned int c = 0; c < send_sm_ranks.size(); c++)
          {
            int i;
            MPI_Waitany(send_sm_req.size(),
                        send_sm_req.data(),
                        &i,
                        MPI_STATUS_IGNORE);

            const Number *__restrict__ data_others_ptr =
              data_others[send_sm_ranks[i]];
            Number *__restrict__ data_this_ptr = data_this;

            for (unsigned int j = send_sm_ptr[i], k = send_sm_offset[i];
                 j < send_sm_ptr[i + 1];
                 j++, k++)
              {
                data_this_ptr[send_sm_indices[j]] += data_others_ptr[k];
              }
          }

        for (unsigned int c = 0; c < send_remote_ranks.size(); c++)
          {
            int i;
            MPI_Waitany(send_remote_req.size(),
                        send_remote_req.data(),
                        &i,
                        MPI_STATUS_IGNORE);

            for (unsigned int j = send_remote_ptr[i];
                 j < send_remote_ptr[i + 1];
                 j++)
              data_this[send_remote_indices[j]] += buffer[j];
          }

        MPI_Waitall(recv_sm_req.size(),
                    recv_sm_req.data(),
                    MPI_STATUSES_IGNORE);

        MPI_Waitall(recv_remote_req.size(),
                    recv_remote_req.data(),
                    MPI_STATUSES_IGNORE);
      }

      void
      compress(Number *data_this, std::vector<Number *> &data_others) const
      {
        compress_start(data_this, data_others);
        compress_finish(data_this, data_others);
      }

    private:
      const MPI_Comm &comm;
      const MPI_Comm &comm_sm;

      unsigned int n_local_elements;

      mutable AlignedVector<Number> buffer;

      std::vector<unsigned int>            recv_remote_ranks;
      std::vector<types::global_dof_index> recv_remote_ptr = {0};
      mutable std::vector<MPI_Request>     recv_remote_req;

      std::vector<unsigned int>        recv_sm_ranks;
      std::vector<unsigned int>        recv_sm_ptr = {0};
      mutable std::vector<MPI_Request> recv_sm_req;
      std::vector<unsigned int>        recv_sm_indices;
      std::vector<unsigned int>        recv_sm_offset;

      std::vector<unsigned int>        send_remote_ranks;
      std::vector<unsigned int>        send_remote_ptr = {0};
      std::vector<unsigned int>        send_remote_indices;
      mutable std::vector<MPI_Request> send_remote_req;

      std::vector<unsigned int>        send_sm_ranks;
      std::vector<unsigned int>        send_sm_ptr = {0};
      std::vector<unsigned int>        send_sm_indices;
      mutable std::vector<MPI_Request> send_sm_req;
      std::vector<unsigned int>        send_sm_offset;
    };

  } // end of namespace SharedMPI
} // end of namespace LinearAlgebra


DEAL_II_NAMESPACE_CLOSE

#endif