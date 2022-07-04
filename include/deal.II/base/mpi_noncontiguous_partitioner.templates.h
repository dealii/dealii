// ---------------------------------------------------------------------
//
// Copyright (C) 2020 - 2022 by the deal.II authors
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

#ifndef dealii_mpi_noncontiguous_partitioner_templates_h
#define dealii_mpi_noncontiguous_partitioner_templates_h

#include <deal.II/base/config.h>

#include <deal.II/base/mpi.h>
#include <deal.II/base/mpi.templates.h>
#include <deal.II/base/mpi_compute_index_owner_internal.h>
#include <deal.II/base/mpi_noncontiguous_partitioner.h>

#include <deal.II/lac/vector_space_vector.h>


DEAL_II_NAMESPACE_OPEN

namespace Utilities
{
  namespace MPI
  {
    template <typename Number>
    void
    NoncontiguousPartitioner::export_to_ghosted_array(
      const ArrayView<const Number> &src,
      const ArrayView<Number> &      dst) const
    {
      // allocate internal memory since needed
      if (requests.size() != send_ranks.size() + recv_ranks.size())
        requests.resize(send_ranks.size() + recv_ranks.size());

      if (this->buffers.size() != send_ptr.back() * sizeof(Number))
        this->buffers.resize(this->temporary_storage_size() * sizeof(Number));

      // perform actual exchange
      this->template export_to_ghosted_array<Number>(
        0,
        src,
        ArrayView<Number>(reinterpret_cast<Number *>(this->buffers.data()),
                          send_ptr.back()),
        dst,
        this->requests);
    }


    template <typename Number>
    void
    NoncontiguousPartitioner::export_to_ghosted_array(
      const unsigned int             communication_channel,
      const ArrayView<const Number> &locally_owned_array,
      const ArrayView<Number> &      temporary_storage,
      const ArrayView<Number> &      ghost_array,
      std::vector<MPI_Request> &     requests) const
    {
      this->template export_to_ghosted_array_start<Number>(
        communication_channel,
        locally_owned_array,
        temporary_storage,
        requests);
      this->template export_to_ghosted_array_finish<Number>(temporary_storage,
                                                            ghost_array,
                                                            requests);
    }



    template <typename Number>
    void
    NoncontiguousPartitioner::export_to_ghosted_array_start(
      const unsigned int             communication_channel,
      const ArrayView<const Number> &src,
      const ArrayView<Number> &      buffers,
      std::vector<MPI_Request> &     requests) const
    {
#ifndef DEAL_II_WITH_MPI
      (void)communication_channel;
      (void)src;
      (void)buffers;
      (void)requests;
      Assert(false, ExcNeedsMPI());
#else
      AssertDimension(requests.size(), recv_ranks.size() + send_ranks.size());

      const auto tag =
        communication_channel +
        internal::Tags::noncontiguous_partitioner_update_ghost_values_start;

      AssertIndexRange(
        tag,
        internal::Tags::noncontiguous_partitioner_update_ghost_values_end + 1);

      // post recv
      AssertIndexRange(recv_ranks.size(), recv_ptr.size());
      for (types::global_dof_index i = 0; i < recv_ranks.size(); ++i)
        {
          const int ierr =
            MPI_Irecv(buffers.data() + recv_ptr[i],
                      recv_ptr[i + 1] - recv_ptr[i],
                      Utilities::MPI::mpi_type_id_for_type<Number>,
                      recv_ranks[i],
                      tag,
                      communicator,
                      &requests[i + send_ranks.size()]);
          AssertThrowMPI(ierr);
        }

      // post send
      AssertIndexRange(send_ranks.size(), send_ptr.size());
      for (types::global_dof_index i = 0, k = 0; i < send_ranks.size(); ++i)
        {
          // collect data to be send
          for (types::global_dof_index j = send_ptr[i]; j < send_ptr[i + 1];
               j++)
            {
              AssertIndexRange(k, send_indices.size());
              buffers[j] = src[send_indices[k]];
              ++k;
            }

          // send data
          Assert((send_ptr[i] < buffers.size()) ||
                   (send_ptr[i] == buffers.size() &&
                    send_ptr[i + 1] == send_ptr[i]),
                 ExcMessage("The input buffer doesn't contain enough entries"));
          const int ierr =
            MPI_Isend(buffers.data() + send_ptr[i],
                      send_ptr[i + 1] - send_ptr[i],
                      Utilities::MPI::mpi_type_id_for_type<Number>,
                      send_ranks[i],
                      tag,
                      communicator,
                      &requests[i]);
          AssertThrowMPI(ierr);
        }
#endif
    }



    template <typename Number>
    void
    NoncontiguousPartitioner::export_to_ghosted_array_finish(
      const ArrayView<const Number> &buffers,
      const ArrayView<Number> &      dst,
      std::vector<MPI_Request> &     requests) const
    {
#ifndef DEAL_II_WITH_MPI
      (void)buffers;
      (void)dst;
      (void)requests;
      Assert(false, ExcNeedsMPI());
#else
      // receive all data packages and copy data from buffers
      for (types::global_dof_index proc = 0; proc < recv_ranks.size(); ++proc)
        {
          int        i;
          MPI_Status status;
          const auto ierr = MPI_Waitany(recv_ranks.size(),
                                        requests.data() + send_ranks.size(),
                                        &i,
                                        &status);
          AssertThrowMPI(ierr);

          AssertIndexRange(i + 1, recv_ptr.size());
          for (types::global_dof_index j = recv_ptr[i], c = 0;
               j < recv_ptr[i + 1];
               j++)
            dst[recv_indices[j]] = buffers[recv_ptr[i] + c++];
        }

      // wait that all data packages have been sent
      const int ierr =
        MPI_Waitall(send_ranks.size(), requests.data(), MPI_STATUSES_IGNORE);
      AssertThrowMPI(ierr);
#endif
    }

  } // namespace MPI
} // namespace Utilities

DEAL_II_NAMESPACE_CLOSE

#endif
