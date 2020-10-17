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

#include <deal.II/base/mpi_compute_index_owner_internal.h>

#include <deal.II/matrix_free/vector_data_exchange.h>

#define DO_COMPRESS true

DEAL_II_NAMESPACE_OPEN

namespace internal
{
  namespace MatrixFreeFunctions
  {
    namespace VectorDataExchange
    {
      namespace internal
      {
        void
        compress(std::vector<unsigned int> &recv_sm_ptr,
                 std::vector<unsigned int> &recv_sm_indices,
                 std::vector<unsigned int> &recv_sm_len)
        {
          std::vector<unsigned int> recv_ptr = {0};
          std::vector<unsigned int> recv_indices;
          std::vector<unsigned int> recv_len;

          for (unsigned int i = 0; i + 1 < recv_sm_ptr.size(); i++)
            {
              if (recv_sm_ptr[i] != recv_sm_ptr[i + 1])
                {
                  recv_indices.push_back(recv_sm_indices[recv_sm_ptr[i]]);
                  recv_len.push_back(1);

                  for (unsigned int j = recv_sm_ptr[i] + 1;
                       j < recv_sm_ptr[i + 1];
                       j++)
                    if (recv_indices.back() + recv_len.back() !=
                        recv_sm_indices[j])
                      {
                        recv_indices.push_back(recv_sm_indices[j]);
                        recv_len.push_back(1);
                      }
                    else
                      recv_len.back()++;
                }
              recv_ptr.push_back(recv_indices.size());
            }

          recv_sm_ptr = recv_ptr;
          recv_sm_ptr.shrink_to_fit();
          recv_sm_indices = recv_indices;
          recv_sm_indices.shrink_to_fit();
          recv_sm_len = recv_len;
          recv_sm_len.shrink_to_fit();
        }
      } // namespace internal


      Full::Full(const IndexSet &is_locally_owned,
                 const IndexSet &is_locally_ghost,
                 const MPI_Comm &communicator,
                 const MPI_Comm &communicator_sm)
      {
        this->comm    = communicator;
        this->comm_sm = communicator_sm;

#ifndef DEAL_II_WITH_MPI
        Assert(false, ExcNeedsMPI());

        (void)is_locally_owned;
        (void)is_locally_ghost;
#else
        this->n_local_elements = is_locally_owned.n_elements();
        this->n_ghost_elements = is_locally_ghost.n_elements();

        this->n_mpi_processes_ = Utilities::MPI::n_mpi_processes(comm);

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

        std::vector<MPI_Request> recv_sm_req;
        std::vector<MPI_Request> send_sm_req;

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

#  if DO_COMPRESS
        internal::compress(recv_sm_ptr, recv_sm_indices, recv_sm_len);
#  endif

#  if DO_COMPRESS
        internal::compress(send_remote_ptr,
                           send_remote_indices,
                           send_remote_len);
        send_remote_offset.clear();
        send_remote_offset.push_back(0);

        for (unsigned int r = 0, c = 0; r < send_remote_ranks.size(); r++)
          {
            for (unsigned int i = send_remote_ptr[r];
                 i < send_remote_ptr[r + 1];
                 i++)
              c += send_remote_len[i];
            send_remote_offset.push_back(c);
          }
#  else
        send_remote_offset = send_remote_ptr;
#  endif

#  if DO_COMPRESS
        internal::compress(send_sm_ptr, send_sm_indices, send_sm_len);
#  endif

#endif
      }

      void
      Full::export_to_ghosted_array_start(
        const unsigned int                          communication_channel,
        const ArrayView<const double> &             locally_owned_array,
        const std::vector<ArrayView<const double>> &shared_arrays,
        const ArrayView<double> &                   ghost_array,
        const ArrayView<double> &                   temporary_storage,
        std::vector<MPI_Request> &                  requests) const
      {
        export_to_ghosted_array_start_impl(communication_channel,
                                           locally_owned_array,
                                           shared_arrays,
                                           ghost_array,
                                           temporary_storage,
                                           requests);
      }

      void
      Full::export_to_ghosted_array_finish(
        const ArrayView<const double> &             locally_owned_array,
        const std::vector<ArrayView<const double>> &shared_arrays,
        const ArrayView<double> &                   ghost_array,
        std::vector<MPI_Request> &                  requests) const
      {
        export_to_ghosted_array_finish_impl(locally_owned_array,
                                            shared_arrays,
                                            ghost_array,
                                            requests);
      }

      void
      Full::import_from_ghosted_array_start(
        const VectorOperation::values               vector_operation,
        const unsigned int                          communication_channel,
        const ArrayView<const double> &             locally_owned_array,
        const std::vector<ArrayView<const double>> &shared_arrays,
        const ArrayView<double> &                   ghost_array,
        const ArrayView<double> &                   temporary_storage,
        std::vector<MPI_Request> &                  requests) const
      {
        import_from_ghosted_array_start_impl(vector_operation,
                                             communication_channel,
                                             locally_owned_array,
                                             shared_arrays,
                                             ghost_array,
                                             temporary_storage,
                                             requests);
      }

      void
      Full::import_from_ghosted_array_finish(
        const VectorOperation::values               vector_operation,
        const ArrayView<double> &                   locally_owned_storage,
        const std::vector<ArrayView<const double>> &shared_arrays,
        const ArrayView<double> &                   ghost_array,
        const ArrayView<const double> &             temporary_storage,
        std::vector<MPI_Request> &                  requests) const
      {
        import_from_ghosted_array_finish_impl(vector_operation,
                                              locally_owned_storage,
                                              shared_arrays,
                                              ghost_array,
                                              temporary_storage,
                                              requests);
      }

      void
      Full::export_to_ghosted_array_start(
        const unsigned int                         communication_channel,
        const ArrayView<const float> &             locally_owned_array,
        const std::vector<ArrayView<const float>> &shared_arrays,
        const ArrayView<float> &                   ghost_array,
        const ArrayView<float> &                   temporary_storage,
        std::vector<MPI_Request> &                 requests) const
      {
        export_to_ghosted_array_start_impl(communication_channel,
                                           locally_owned_array,
                                           shared_arrays,
                                           ghost_array,
                                           temporary_storage,
                                           requests);
      }

      void
      Full::export_to_ghosted_array_finish(
        const ArrayView<const float> &             locally_owned_array,
        const std::vector<ArrayView<const float>> &shared_arrays,
        const ArrayView<float> &                   ghost_array,
        std::vector<MPI_Request> &                 requests) const
      {
        export_to_ghosted_array_finish_impl(locally_owned_array,
                                            shared_arrays,
                                            ghost_array,
                                            requests);
      }

      void
      Full::import_from_ghosted_array_start(
        const VectorOperation::values              vector_operation,
        const unsigned int                         communication_channel,
        const ArrayView<const float> &             locally_owned_array,
        const std::vector<ArrayView<const float>> &shared_arrays,
        const ArrayView<float> &                   ghost_array,
        const ArrayView<float> &                   temporary_storage,
        std::vector<MPI_Request> &                 requests) const
      {
        import_from_ghosted_array_start_impl(vector_operation,
                                             communication_channel,
                                             locally_owned_array,
                                             shared_arrays,
                                             ghost_array,
                                             temporary_storage,
                                             requests);
      }

      void
      Full::import_from_ghosted_array_finish(
        const VectorOperation::values              vector_operation,
        const ArrayView<float> &                   locally_owned_storage,
        const std::vector<ArrayView<const float>> &shared_arrays,
        const ArrayView<float> &                   ghost_array,
        const ArrayView<const float> &             temporary_storage,
        std::vector<MPI_Request> &                 requests) const
      {
        import_from_ghosted_array_finish_impl(vector_operation,
                                              locally_owned_storage,
                                              shared_arrays,
                                              ghost_array,
                                              temporary_storage,
                                              requests);
      }

      template <typename Number>
      void
      Full::export_to_ghosted_array_start_impl(
        const unsigned int                          communication_channel,
        const ArrayView<const Number> &             data_this,
        const std::vector<ArrayView<const Number>> &data_others,
        const ArrayView<Number> &                   buffer,
        const ArrayView<Number> &                   temporary_storage,
        std::vector<MPI_Request> &                  requests) const
      {
        (void)data_this;

#ifndef DEAL_II_WITH_MPI
        Assert(false, ExcNeedsMPI());

        (void)temporary_storage;
        (void)communication_channel;
        (void)data_others;
        (void)buffer;
        (void)requests;
#else
        (void)data_others;

        requests.resize(send_sm_ranks.size() + recv_sm_ranks.size() +
                        recv_remote_ranks.size() + send_remote_ranks.size());

        int dummy;
        // receive a signal that relevant sm neighbors are ready
        for (unsigned int i = 0; i < recv_sm_ranks.size(); i++)
          MPI_Irecv(&dummy,
                    0,
                    MPI_INT,
                    recv_sm_ranks[i],
                    communication_channel + 2,
                    comm_sm,
                    requests.data() + send_sm_ranks.size() + i);

        // signal to all relevant sm neighbors that this process is ready
        for (unsigned int i = 0; i < send_sm_ranks.size(); i++)
          MPI_Isend(&dummy,
                    0,
                    MPI_INT,
                    send_sm_ranks[i],
                    communication_channel + 2,
                    comm_sm,
                    requests.data() + i);

        // receive data from remote processes
        for (unsigned int i = 0; i < recv_remote_ranks.size(); i++)
          MPI_Irecv(buffer.data() + recv_remote_ptr[i],
                    recv_remote_ptr[i + 1] - recv_remote_ptr[i],
                    Utilities::MPI::internal::mpi_type_id(data_this.data()),
                    recv_remote_ranks[i],
                    communication_channel + 3,
                    comm,
                    requests.data() + send_sm_ranks.size() +
                      recv_sm_ranks.size() + i);

          // send data to remote processes
#  if DO_COMPRESS
        for (unsigned int i = 0, k = 0; i < send_remote_ranks.size(); i++)
          {
            for (unsigned int j = send_remote_ptr[i];
                 j < send_remote_ptr[i + 1];
                 j++)
              for (unsigned int l = 0; l < send_remote_len[j]; l++, k++)
                temporary_storage[k] = data_this[send_remote_indices[j] + l];
#  else
        for (unsigned int i = 0; i < send_remote_ranks.size(); i++)
          {
            for (unsigned int j = send_remote_ptr[i];
                 j < send_remote_ptr[i + 1];
                 j++)
              temporary_storage[j] = data_this[send_remote_indices[j]];
#  endif

            // send data away
            MPI_Isend(temporary_storage.data() + send_remote_offset[i],
                      send_remote_offset[i + 1] - send_remote_offset[i],
                      Utilities::MPI::internal::mpi_type_id(data_this.data()),
                      send_remote_ranks[i],
                      communication_channel + 3,
                      comm,
                      requests.data() + send_sm_ranks.size() +
                        recv_sm_ranks.size() + recv_remote_ranks.size() + i);
          }
#endif
      }

      template <typename Number>
      void
      Full::export_to_ghosted_array_finish_impl(
        const ArrayView<const Number> &             data_this,
        const std::vector<ArrayView<const Number>> &data_others,
        const ArrayView<Number> &                   ghost_array,
        std::vector<MPI_Request> &                  requests) const
      {
        (void)data_this;

#ifndef DEAL_II_WITH_MPI
        Assert(false, ExcNeedsMPI());

        (void)ghost_array;
        (void)data_others;
        (void)requests;
#else
        AssertDimension(requests.size(),
                        send_sm_ranks.size() + recv_sm_ranks.size() +
                          recv_remote_ranks.size() + send_remote_ranks.size());

        for (unsigned int c = 0; c < recv_sm_ranks.size(); c++)
          {
            int i;
            MPI_Waitany(recv_sm_ranks.size(),
                        requests.data() + send_sm_ranks.size(),
                        &i,
                        MPI_STATUS_IGNORE);

            const Number *__restrict__ data_others_ptr =
              data_others[recv_sm_ranks[i]].data();
            Number *__restrict__ data_this_ptr = ghost_array.data();

#  if DO_COMPRESS
            for (unsigned int j = recv_sm_ptr[i], k = recv_sm_offset[i];
                 j < recv_sm_ptr[i + 1];
                 j++)
              for (unsigned int l = 0; l < recv_sm_len[j]; l++, k++)
                data_this_ptr[k] =
                  data_others_ptr[recv_sm_indices[j] + l]; // TODO!!!
#  else
            for (unsigned int j = recv_sm_ptr[i], k = recv_sm_offset[i];
                 j < recv_sm_ptr[i + 1];
                 j++, k++)
              data_this_ptr[k] = data_others_ptr[recv_sm_indices[j]]; // TODO!!!
#  endif
          }

        MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
#endif
      }

      template <typename Number>
      void
      Full::import_from_ghosted_array_start_impl(
        const VectorOperation::values               operation,
        const unsigned int                          communication_channel,
        const ArrayView<const Number> &             data_this,
        const std::vector<ArrayView<const Number>> &data_others,
        const ArrayView<Number> &                   buffer,
        const ArrayView<Number> &                   temporary_storage,
        std::vector<MPI_Request> &                  requests) const
      {
        (void)data_this;

#ifndef DEAL_II_WITH_MPI
        Assert(false, ExcNeedsMPI());

        (void)operation;
        (void)communication_channel;
        (void)data_others;
        (void)buffer;
        (void)temporary_storage;
        (void)requests;
#else
        (void)data_others;
        (void)operation;

        Assert(operation == dealii::VectorOperation::add, ExcNotImplemented());

        requests.resize(recv_sm_ranks.size() + send_sm_ranks.size() +
                        recv_remote_ranks.size() + send_remote_ranks.size());

        int dummy;
        for (unsigned int i = 0; i < recv_sm_ranks.size(); i++)
          MPI_Isend(&dummy,
                    0,
                    MPI_INT,
                    recv_sm_ranks[i],
                    communication_channel + 1,
                    comm_sm,
                    requests.data() + i);

        for (unsigned int i = 0; i < send_sm_ranks.size(); i++)
          MPI_Irecv(&dummy,
                    0,
                    MPI_INT,
                    send_sm_ranks[i],
                    communication_channel + 1,
                    comm_sm,
                    requests.data() + recv_sm_ranks.size() + i);

        for (unsigned int i = 0; i < recv_remote_ranks.size(); i++)
          MPI_Isend(buffer.data() + recv_remote_ptr[i],
                    recv_remote_ptr[i + 1] - recv_remote_ptr[i],
                    Utilities::MPI::internal::mpi_type_id(buffer.data()),
                    recv_remote_ranks[i],
                    communication_channel + 0,
                    comm,
                    requests.data() + recv_sm_ranks.size() +
                      send_sm_ranks.size() + i);

        for (unsigned int i = 0; i < send_remote_ranks.size(); i++)
          MPI_Irecv(temporary_storage.data() + send_remote_offset[i],
                    send_remote_offset[i + 1] - send_remote_offset[i],
                    Utilities::MPI::internal::mpi_type_id(
                      temporary_storage.data()),
                    send_remote_ranks[i],
                    communication_channel + 0,
                    comm,
                    requests.data() + recv_sm_ranks.size() +
                      send_sm_ranks.size() + recv_remote_ranks.size() + i);
#endif
      }

      template <typename Number>
      void
      Full::import_from_ghosted_array_finish_impl(
        const VectorOperation::values               operation,
        const ArrayView<Number> &                   data_this,
        const std::vector<ArrayView<const Number>> &data_others,
        const ArrayView<Number> &                   buffer,
        const ArrayView<const Number> &             temporary_storage,
        std::vector<MPI_Request> &                  requests) const
      {
        (void)temporary_storage;

#ifndef DEAL_II_WITH_MPI
        Assert(false, ExcNeedsMPI());

        (void)operation;
        (void)data_this;
        (void)data_others;
        (void)buffer;
        (void)requests;
#else
        (void)operation;

        Assert(operation == dealii::VectorOperation::add, ExcNotImplemented());

        AssertDimension(requests.size(),
                        recv_sm_ranks.size() + send_sm_ranks.size() +
                          recv_remote_ranks.size() + send_remote_ranks.size());

        const auto split = [&](const unsigned int i) {
          AssertIndexRange(i,
                           (send_sm_ranks.size() + recv_remote_ranks.size() +
                            send_remote_ranks.size()));

          if (i < send_sm_ranks.size())
            return std::pair<unsigned int, unsigned int>{0, i};
          else if (i < (send_sm_ranks.size() + recv_remote_ranks.size()))
            return std::pair<unsigned int, unsigned int>{
              2, i - send_sm_ranks.size()};
          else
            return std::pair<unsigned int, unsigned int>{
              1, i - send_sm_ranks.size() - recv_remote_ranks.size()};
        };

        for (unsigned int c = 0;
             c < send_sm_ranks.size() + send_remote_ranks.size() +
                   recv_remote_ranks.size();
             c++)
          {
            int i;
            MPI_Waitany(send_sm_ranks.size() + send_remote_ranks.size() +
                          recv_remote_ranks.size(),
                        requests.data() + recv_sm_ranks.size(),
                        &i,
                        MPI_STATUS_IGNORE);

            const auto &s = split(i);
            i             = s.second;

            if (s.first == 0)
              {
                Number *__restrict__ data_others_ptr =
                  const_cast<Number *>(data_others[send_sm_ranks[i]].data());
                Number *__restrict__ data_this_ptr = data_this.data();

#  if DO_COMPRESS
                for (unsigned int j = send_sm_ptr[i], k = send_sm_offset[i];
                     j < send_sm_ptr[i + 1];
                     j++)
                  {
                    for (unsigned int l = 0; l < send_sm_len[j]; l++, k++)
                      {
                        data_this_ptr[send_sm_indices[j] + l] +=
                          data_others_ptr[k];
                        data_others_ptr[k] = 0.0;
                      }
                  }
#  else
                for (unsigned int j = send_sm_ptr[i], k = send_sm_offset[i];
                     j < send_sm_ptr[i + 1];
                     j++, k++)
                  {
                    data_this_ptr[send_sm_indices[j]] += data_others_ptr[k];
                    data_others_ptr[k] = 0.0;
                  }
#  endif
              }
            else if (s.first == 1)
              {
#  if DO_COMPRESS
                for (unsigned int j = send_remote_ptr[i],
                                  k = send_remote_offset[i];
                     j < send_remote_ptr[i + 1];
                     j++)
                  for (unsigned int l = 0; l < send_remote_len[j]; l++)
                    data_this[send_remote_indices[j] + l] += buffer[k++];
#  else
                for (unsigned int j = send_remote_ptr[i];
                     j < send_remote_ptr[i + 1];
                     j++)
                  data_this[send_remote_indices[j]] += buffer[j];
#  endif
              }
            else /*if (s.first == 2)*/
              {
                std::memset(data_this.data() + recv_remote_ptr[i] +
                              n_local_elements,
                            0,
                            (recv_remote_ptr[i + 1] - recv_remote_ptr[i]) *
                              sizeof(Number));
              }
          }

        MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
#endif
      }

      unsigned int
      Full::local_size() const
      {
        return n_local_elements;
      }

      unsigned int
      Full::n_ghost_indices() const
      {
        return n_ghost_elements;
      }

      unsigned int
      Full::n_import_indices() const
      {
        return send_remote_offset.back();
      }

      void
      Full::reset_ghost_values(const ArrayView<double> &ghost_array) const
      {
        (void)ghost_array;
        // nothing to do
      }

      void
      Full::reset_ghost_values(const ArrayView<float> &ghost_array) const
      {
        (void)ghost_array;
        // nothing to do
      }

    } // namespace VectorDataExchange
  }   // namespace MatrixFreeFunctions
} // namespace internal


DEAL_II_NAMESPACE_CLOSE
