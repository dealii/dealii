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

#include <deal.II/lac/la_sm_partitioner.h>

#define DO_COMPRESS true

DEAL_II_NAMESPACE_OPEN

namespace LinearAlgebra
{
  namespace SharedMPI
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

        //        for(auto i : recv_sm_ptr)
        //            std::cout << i << " ";
        //        std::cout << std::endl;
        //
        //        for(auto i : recv_sm_indices)
        //            std::cout << i << " ";
        //        std::cout << std::endl;

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

        //        for(auto i : recv_sm_ptr)
        //            std::cout << i << " ";
        //        std::cout << std::endl;
        //
        //        for(auto i : recv_sm_indices)
        //            std::cout << i << " ";
        //        std::cout << std::endl;
        //
        //        for(auto i : recv_sm_len)
        //            std::cout << i << " ";
        //        std::cout << std::endl;
      }

      template <typename Number = double>
      std::vector<unsigned int>
      create_sm_view(const MPI_Comm &           comm,
                     const unsigned int         len,
                     std::vector<unsigned int> &recv_sm_ranks,
                     std::vector<unsigned int> &recv_sm_ptr,
                     std::vector<unsigned int> &recv_sm_indices,
                     std::vector<unsigned int> &recv_sm_offset)
      {
        std::vector<unsigned int> offsets(
          Utilities::MPI::n_mpi_processes(comm) + 1);

        const unsigned int my_rank = Utilities::MPI::this_mpi_process(comm);

        unsigned int len_ = ((len * sizeof(Number) + 64 - 1) / sizeof(Number));

        MPI_Allgather(&len_,
                      1,
                      Utilities::MPI::internal::mpi_type_id(offsets.data()),
                      offsets.data() + 1,
                      1,
                      Utilities::MPI::internal::mpi_type_id(offsets.data()),
                      comm);

        for (unsigned int i = 1; i < offsets.size(); i++)
          offsets[i] += offsets[i - 1];

        //        for(auto i : offsets)
        //            std::cout << i << " ";
        //        std::cout << std::endl;

        std::vector<unsigned int> recv_sm_indices_temp(
          len, numbers::invalid_unsigned_int);

        for (unsigned int i = 0; i < recv_sm_ranks.size(); i++)
          for (unsigned int j = recv_sm_ptr[i], c = 0; j < recv_sm_ptr[i + 1];
               j++, c++)
            recv_sm_indices_temp[recv_sm_offset[i] + c] =
              recv_sm_indices[j] + offsets[recv_sm_ranks[i]];

        for (unsigned int i = 0; i < len; i++)
          if (recv_sm_indices_temp[i] == numbers::invalid_unsigned_int)
            recv_sm_indices_temp[i] = i + offsets[my_rank];


        //        for(auto i : recv_sm_indices)
        //            std::cout << i << " ";
        //        std::cout << std::endl;
        //        for(auto i : recv_sm_indices_temp)
        //            std::cout << i << " ";
        //        std::cout << std::endl;

        return recv_sm_indices_temp;
      }

    } // namespace internal

    Partitioner::Partitioner(const MPI_Comm &comm,
                             const MPI_Comm &comm_sm,
                             const IndexSet &is_locally_owned,
                             const IndexSet &is_locally_ghost)
      : comm(comm)
      , comm_sm(comm_sm)
    {
      reinit(is_locally_owned, is_locally_ghost);
    }

    const MPI_Comm &
    Partitioner::get_mpi_communicator() const
    {
      return comm;
    }

    const MPI_Comm &
    Partitioner::get_sm_mpi_communicator() const
    {
      return comm_sm;
    }

    void
    Partitioner::reinit(const IndexSet &is_locally_owned,
                        const IndexSet &is_locally_ghost,
                        const MPI_Comm &communicator)
    {
      Assert(false, ExcNotImplemented());
      (void)is_locally_owned;
      (void)is_locally_ghost;
      (void)communicator;
    }

    void
    Partitioner::reinit(const IndexSet &is_locally_owned,
                        const IndexSet &is_locally_ghost)
    {
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
                recv_remote_ptr.push_back(recv_remote_ptr.back() +
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

      //      const auto a =
      //        Utilities::MPI::sum(static_cast<double>(this->memory_consumption()),
      //                            this->comm);

      this->sm_view =
        internal::create_sm_view(comm_sm,
                                 this->local_size() + this->n_ghost_indices(),
                                 recv_sm_ranks,
                                 recv_sm_ptr,
                                 recv_sm_indices,
                                 recv_sm_offset);

#if DO_COMPRESS
      internal::compress(recv_sm_ptr, recv_sm_indices, recv_sm_len);
#endif

#if DO_COMPRESS
      internal::compress(send_remote_ptr, send_remote_indices, send_remote_len);
      send_remote_offset.clear();
      send_remote_offset.push_back(0);

      for (unsigned int r = 0, c = 0; r < send_remote_ranks.size(); r++)
        {
          for (unsigned int i = send_remote_ptr[r]; i < send_remote_ptr[r + 1];
               i++)
            c += send_remote_len[i];
          send_remote_offset.push_back(c);
        }
#else
      send_remote_offset = send_remote_ptr;
#endif

#if DO_COMPRESS
      internal::compress(send_sm_ptr, send_sm_indices, send_sm_len);
#endif

      //      const auto b =
      //        Utilities::MPI::sum(static_cast<double>(this->memory_consumption()),
      //                            this->comm);
      //
      //      std::cout << "memory_consumption " << a << " " << b << std::endl;
    }

    void
    Partitioner::export_to_ghosted_array_start(
      const unsigned int             communication_channel,
      double *const                  data_this,
      const std::vector<double *> &  data_others,
      dealii::AlignedVector<double> &buffer) const
    {
      export_to_ghosted_array_start_impl(communication_channel,
                                         data_this,
                                         data_others,
                                         buffer);
    }

    void
    Partitioner::export_to_ghosted_array_finish(
      double *const                data_this,
      const std::vector<double *> &data_others) const
    {
      export_to_ghosted_array_finish_impl(data_this, data_others);
    }

    void
    Partitioner::import_from_ghosted_array_start(
      const VectorOperation::values  operation,
      const unsigned int             communication_channel,
      double *const                  data_this,
      const std::vector<double *> &  data_others,
      dealii::AlignedVector<double> &buffer) const
    {
      import_from_ghosted_array_start_impl(
        operation, communication_channel, data_this, data_others, buffer);
    }

    void
    Partitioner::import_from_ghosted_array_finish(
      const VectorOperation::values        operation,
      double *const                        data_this,
      const std::vector<double *> &        data_others,
      const dealii::AlignedVector<double> &buffer) const
    {
      import_from_ghosted_array_finish_impl(operation,
                                            data_this,
                                            data_others,
                                            buffer);
    }

    void
    Partitioner::export_to_ghosted_array_start(
      const unsigned int            communication_channel,
      float *const                  data_this,
      const std::vector<float *> &  data_others,
      dealii::AlignedVector<float> &buffer) const
    {
      export_to_ghosted_array_start_impl(communication_channel,
                                         data_this,
                                         data_others,
                                         buffer);
    }

    void
    Partitioner::export_to_ghosted_array_finish(
      float *const                data_this,
      const std::vector<float *> &data_others) const
    {
      export_to_ghosted_array_finish_impl(data_this, data_others);
    }

    void
    Partitioner::import_from_ghosted_array_start(
      const VectorOperation::values operation,
      const unsigned int            communication_channel,
      float *const                  data_this,
      const std::vector<float *> &  data_others,
      dealii::AlignedVector<float> &buffer) const
    {
      import_from_ghosted_array_start_impl(
        operation, communication_channel, data_this, data_others, buffer);
    }

    void
    Partitioner::import_from_ghosted_array_finish(
      const VectorOperation::values       operation,
      float *const                        data_this,
      const std::vector<float *> &        data_others,
      const dealii::AlignedVector<float> &buffer) const
    {
      import_from_ghosted_array_finish_impl(operation,
                                            data_this,
                                            data_others,
                                            buffer);
    }

    template <typename Number>
    void
    Partitioner::export_to_ghosted_array_start_impl(
      const unsigned int             communication_channel,
      Number *const                  data_this,
      const std::vector<Number *> &  data_others,
      dealii::AlignedVector<Number> &buffer) const
    {
      (void)data_others;

      if (send_remote_offset.back() != buffer.size())
        buffer.resize(send_remote_offset.back());

      int dummy;
      for (unsigned int i = 0; i < recv_sm_ranks.size(); i++)
        MPI_Irecv(&dummy,
                  0,
                  MPI_INT,
                  recv_sm_ranks[i],
                  communication_channel + 2,
                  comm_sm,
                  recv_sm_req.data() + i);

      for (unsigned int i = 0; i < send_sm_ranks.size(); i++)
        MPI_Isend(&dummy,
                  0,
                  MPI_INT,
                  send_sm_ranks[i],
                  communication_channel + 2,
                  comm_sm,
                  send_sm_req.data() + i);

      for (unsigned int i = 0; i < recv_remote_ranks.size(); i++)
        MPI_Irecv(data_this + recv_remote_ptr[i] + n_local_elements,
                  recv_remote_ptr[i + 1] - recv_remote_ptr[i],
                  Utilities::MPI::internal::mpi_type_id(data_this),
                  recv_remote_ranks[i],
                  communication_channel + 3,
                  comm,
                  recv_remote_req.data() + i);

#if DO_COMPRESS
      for (unsigned int i = 0, k = 0; i < send_remote_ranks.size(); i++)
        {
          for (unsigned int j = send_remote_ptr[i]; j < send_remote_ptr[i + 1];
               j++)
            for (unsigned int l = 0; l < send_remote_len[j]; l++, k++)
              buffer[k] = data_this[send_remote_indices[j] + l];
#else
      for (unsigned int i = 0; i < send_remote_ranks.size(); i++)
        {
          for (unsigned int j = send_remote_ptr[i]; j < send_remote_ptr[i + 1];
               j++)
            buffer[j] = data_this[send_remote_indices[j]];
#endif

          MPI_Isend(buffer.data() + send_remote_offset[i],
                    send_remote_offset[i + 1] - send_remote_offset[i],
                    Utilities::MPI::internal::mpi_type_id(data_this),
                    send_remote_ranks[i],
                    communication_channel + 3,
                    comm,
                    send_remote_req.data() + i);
        }
    }

    template <typename Number>
    void
    Partitioner::export_to_ghosted_array_finish_impl(
      Number *const                data_this,
      const std::vector<Number *> &data_others) const
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

#if DO_COMPRESS
          for (unsigned int j = recv_sm_ptr[i], k = recv_sm_offset[i];
               j < recv_sm_ptr[i + 1];
               j++)
            for (unsigned int l = 0; l < recv_sm_len[j]; l++, k++)
              data_this_ptr[k] = data_others_ptr[recv_sm_indices[j] + l];
#else
          for (unsigned int j = recv_sm_ptr[i], k = recv_sm_offset[i];
               j < recv_sm_ptr[i + 1];
               j++, k++)
            data_this_ptr[k] = data_others_ptr[recv_sm_indices[j]];
#endif
        }

      MPI_Waitall(send_sm_req.size(), send_sm_req.data(), MPI_STATUSES_IGNORE);
      MPI_Waitall(send_remote_req.size(),
                  send_remote_req.data(),
                  MPI_STATUSES_IGNORE);
      MPI_Waitall(recv_remote_req.size(),
                  recv_remote_req.data(),
                  MPI_STATUSES_IGNORE);
    }

    template <typename Number>
    void
    Partitioner::import_from_ghosted_array_start_impl(
      const VectorOperation::values  operation,
      const unsigned int             communication_channel,
      Number *const                  data_this,
      const std::vector<Number *> &  data_others,
      dealii::AlignedVector<Number> &buffer) const
    {
      (void)data_others;
      (void)operation;

      Assert(operation == dealii::VectorOperation::add, ExcNotImplemented());

      if (send_remote_offset.back() != buffer.size())
        buffer.resize(send_remote_offset.back());

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
        MPI_Irecv(buffer.data() + send_remote_offset[i],
                  send_remote_offset[i + 1] - send_remote_offset[i],
                  Utilities::MPI::internal::mpi_type_id(data_this),
                  send_remote_ranks[i],
                  communication_channel + 0,
                  comm,
                  send_remote_req.data() + i);
    }

    template <typename Number>
    void
    Partitioner::import_from_ghosted_array_finish_impl(
      const VectorOperation::values        operation,
      Number *const                        data_this,
      const std::vector<Number *> &        data_others,
      const dealii::AlignedVector<Number> &buffer) const
    {
      (void)operation;

      Assert(operation == dealii::VectorOperation::add, ExcNotImplemented());

      for (unsigned int c = 0; c < send_sm_ranks.size(); c++)
        {
          int i;
          MPI_Waitany(send_sm_req.size(),
                      send_sm_req.data(),
                      &i,
                      MPI_STATUS_IGNORE);

          Number *__restrict__ data_others_ptr = data_others[send_sm_ranks[i]];
          Number *__restrict__ data_this_ptr   = data_this;

#if DO_COMPRESS
          for (unsigned int j = send_sm_ptr[i], k = send_sm_offset[i];
               j < send_sm_ptr[i + 1];
               j++)
            {
              for (unsigned int l = 0; l < send_sm_len[j]; l++, k++)
                {
                  data_this_ptr[send_sm_indices[j] + l] += data_others_ptr[k];
                  data_others_ptr[k] = 0.0;
                }
            }
#else
          for (unsigned int j = send_sm_ptr[i], k = send_sm_offset[i];
               j < send_sm_ptr[i + 1];
               j++, k++)
            {
              data_this_ptr[send_sm_indices[j]] += data_others_ptr[k];
              data_others_ptr[k] = 0.0;
            }
#endif
        }

      for (unsigned int c = 0; c < send_remote_ranks.size(); c++)
        {
          int i;
          MPI_Waitany(send_remote_req.size(),
                      send_remote_req.data(),
                      &i,
                      MPI_STATUS_IGNORE);

#if DO_COMPRESS
          for (unsigned int j = send_remote_ptr[i], k = send_remote_offset[i];
               j < send_remote_ptr[i + 1];
               j++)
            for (unsigned int l = 0; l < send_remote_len[j]; l++)
              data_this[send_remote_indices[j] + l] += buffer[k++];
#else
          for (unsigned int j = send_remote_ptr[i]; j < send_remote_ptr[i + 1];
               j++)
            data_this[send_remote_indices[j]] += buffer[j];
#endif
        }

      MPI_Waitall(recv_sm_req.size(), recv_sm_req.data(), MPI_STATUSES_IGNORE);

      MPI_Waitall(recv_remote_req.size(),
                  recv_remote_req.data(),
                  MPI_STATUSES_IGNORE);
    }

    std::size_t
    Partitioner::local_size() const
    {
      return n_local_elements;
    }

    std::size_t
    Partitioner::n_ghost_indices() const
    {
      return n_ghost_elements;
    }

    std::size_t
    Partitioner::n_mpi_processes() const
    {
      return n_mpi_processes_;
    }

    std::size_t
    Partitioner::memory_consumption() const
    {
      return 0                                                            //
             + MemoryConsumption::memory_consumption(recv_remote_ranks)   //
             + MemoryConsumption::memory_consumption(recv_remote_ptr)     //
             + MemoryConsumption::memory_consumption(recv_sm_ranks)       //
             + MemoryConsumption::memory_consumption(recv_sm_ptr)         //
             + MemoryConsumption::memory_consumption(recv_sm_indices)     //
             + MemoryConsumption::memory_consumption(recv_sm_len)         //
             + MemoryConsumption::memory_consumption(recv_sm_offset)      //
             + MemoryConsumption::memory_consumption(send_remote_ptr)     //
             + MemoryConsumption::memory_consumption(send_remote_indices) //
             + MemoryConsumption::memory_consumption(send_remote_len)     //
             + MemoryConsumption::memory_consumption(send_remote_offset)  //
             + MemoryConsumption::memory_consumption(send_sm_ranks)       //
             + MemoryConsumption::memory_consumption(send_sm_ptr)         //
             + MemoryConsumption::memory_consumption(send_sm_indices)     //
             + MemoryConsumption::memory_consumption(send_sm_len)         //
             + MemoryConsumption::memory_consumption(send_sm_offset)      //
        ;
    }

    std::vector<unsigned int>
    Partitioner::get_sm_view() const
    {
      return sm_view;
    }



  } // end of namespace SharedMPI
} // end of namespace LinearAlgebra

#include "la_sm_partitioner.inst"


DEAL_II_NAMESPACE_CLOSE
