// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_mpi_noncontiguous_partitioner_templates_h
#define dealii_mpi_noncontiguous_partitioner_templates_h

#include <deal.II/base/config.h>

#include <deal.II/base/mpi.h>
#include <deal.II/base/mpi.templates.h>
#include <deal.II/base/mpi_noncontiguous_partitioner.h>
#include <deal.II/base/mpi_tags.h>



DEAL_II_NAMESPACE_OPEN

namespace Utilities
{
  namespace MPI
  {
    template <typename Number, unsigned int n_components_templated>
    void
    NoncontiguousPartitioner::export_to_ghosted_array(
      const ArrayView<const Number> &src,
      const ArrayView<Number>       &dst,
      const unsigned int             n_components) const
    {
      Assert((n_components_templated != 0) != (n_components != 0),
             ExcNotImplemented());

      const unsigned int n_components_to_be_used =
        (n_components_templated != 0) ? n_components_templated : n_components;

      // allocate internal memory since needed
      if (requests.size() != send_ranks.size() + recv_ranks.size())
        requests.resize(send_ranks.size() + recv_ranks.size());

      if (this->buffers.size() !=
          send_ptr.back() * sizeof(Number) * n_components_to_be_used)
        this->buffers.resize(this->temporary_storage_size() * sizeof(Number) *
                             n_components_to_be_used);

      // perform actual exchange
      this->template export_to_ghosted_array<Number, n_components_templated>(
        0,
        src,
        ArrayView<Number>(reinterpret_cast<Number *>(this->buffers.data()),
                          send_ptr.back() * n_components_to_be_used),
        dst,
        this->requests,
        n_components);
    }


    template <typename Number, unsigned int n_components_templated>
    void
    NoncontiguousPartitioner::export_to_ghosted_array(
      const unsigned int             communication_channel,
      const ArrayView<const Number> &locally_owned_array,
      const ArrayView<Number>       &temporary_storage,
      const ArrayView<Number>       &ghost_array,
      std::vector<MPI_Request>      &requests,
      const unsigned int             n_components) const
    {
      this->template export_to_ghosted_array_start<Number,
                                                   n_components_templated>(
        communication_channel,
        locally_owned_array,
        temporary_storage,
        requests,
        n_components);
      this->template export_to_ghosted_array_finish<Number,
                                                    n_components_templated>(
        temporary_storage, ghost_array, requests, n_components);
    }



    template <typename Number, unsigned int n_components_templated>
    void
    NoncontiguousPartitioner::export_to_ghosted_array_start(
      const unsigned int             communication_channel,
      const ArrayView<const Number> &src,
      const ArrayView<Number>       &buffers,
      std::vector<MPI_Request>      &requests,
      const unsigned int             n_components) const
    {
#ifndef DEAL_II_WITH_MPI
      (void)communication_channel;
      (void)src;
      (void)buffers;
      (void)requests;
      (void)n_components;
      Assert(false, ExcNeedsMPI());
#else
      AssertDimension(requests.size(), recv_ranks.size() + send_ranks.size());

      Assert((n_components_templated != 0) != (n_components != 0),
             ExcNotImplemented());

      const unsigned int n_components_to_be_used =
        (n_components_templated != 0) ? n_components_templated : n_components;

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
            MPI_Irecv(buffers.data() + recv_ptr[i] * n_components_to_be_used,
                      (recv_ptr[i + 1] - recv_ptr[i]) * n_components_to_be_used,
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
               ++j, ++k)
            {
              AssertIndexRange(k, send_indices.size());
              for (unsigned int comp = 0; comp < n_components_to_be_used;
                   ++comp)
                buffers[j * n_components_to_be_used + comp] =
                  src[send_indices[k] * n_components_to_be_used + comp];
            }

          // send data
          Assert((send_ptr[i] * n_components_to_be_used < buffers.size()) ||
                   (send_ptr[i] * n_components_to_be_used == buffers.size() &&
                    send_ptr[i + 1] * n_components_to_be_used ==
                      send_ptr[i] * n_components_to_be_used),
                 ExcMessage("The input buffer doesn't contain enough entries"));
          const int ierr =
            MPI_Isend(buffers.data() + send_ptr[i] * n_components_to_be_used,
                      (send_ptr[i + 1] - send_ptr[i]) * n_components_to_be_used,
                      Utilities::MPI::mpi_type_id_for_type<Number>,
                      send_ranks[i],
                      tag,
                      communicator,
                      &requests[i]);
          AssertThrowMPI(ierr);
        }
#endif
    }



    template <typename Number, unsigned int n_components_templated>
    void
    NoncontiguousPartitioner::export_to_ghosted_array_finish(
      const ArrayView<const Number> &buffers,
      const ArrayView<Number>       &dst,
      std::vector<MPI_Request>      &requests,
      const unsigned int             n_components) const
    {
#ifndef DEAL_II_WITH_MPI
      (void)buffers;
      (void)dst;
      (void)requests;
      (void)n_components;
      Assert(false, ExcNeedsMPI());
#else

      Assert((n_components_templated != 0) != (n_components != 0),
             ExcNotImplemented());

      const unsigned int n_components_to_be_used =
        (n_components_templated != 0) ? n_components_templated : n_components;

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
          for (types::global_dof_index j = recv_ptr[i]; j < recv_ptr[i + 1];
               ++j)
            for (unsigned int comp = 0; comp < n_components_to_be_used; ++comp)
              {
                const auto &value = buffers[j * n_components_to_be_used + comp];

                if (recv_indices_duplicates_ptr.empty())
                  dst[recv_indices[j] * n_components_to_be_used + comp] = value;
                else
                  for (auto k = recv_indices_duplicates_ptr[recv_indices[j]];
                       k < recv_indices_duplicates_ptr[recv_indices[j] + 1];
                       ++k)
                    dst[recv_indices_duplicates[k] * n_components_to_be_used +
                        comp] = value;
              }
        }

      // wait that all data packages have been sent
      const int ierr =
        MPI_Waitall(send_ranks.size(), requests.data(), MPI_STATUSES_IGNORE);
      AssertThrowMPI(ierr);
#endif
    }



    template <typename Number>
    void
    NoncontiguousPartitioner::import_from_ghosted_array(
      const VectorOperation::values vector_operation,
      const ArrayView<Number>      &src,
      const ArrayView<Number>      &dst) const
    {
      // allocate internal memory since needed
      if (requests.size() != send_ranks.size() + recv_ranks.size())
        requests.resize(send_ranks.size() + recv_ranks.size());

      if (this->buffers.size() != send_ptr.back() * sizeof(Number))
        this->buffers.resize(this->temporary_storage_size() * sizeof(Number));

      // perform actual exchange
      this->template import_from_ghosted_array<Number>(
        vector_operation,
        0,
        src,
        ArrayView<Number>(reinterpret_cast<Number *>(this->buffers.data()),
                          send_ptr.back()),
        dst,
        this->requests);
    }



    template <typename Number>
    void
    NoncontiguousPartitioner::import_from_ghosted_array(
      const VectorOperation::values vector_operation,
      const unsigned int            communication_channel,
      const ArrayView<Number>      &ghost_array,
      const ArrayView<Number>      &temporary_storage,
      const ArrayView<Number>      &locally_owned_array,
      std::vector<MPI_Request>     &requests) const
    {
      this->template import_from_ghosted_array_start<Number>(
        vector_operation,
        communication_channel,
        ghost_array,
        temporary_storage,
        requests);
      this->template import_from_ghosted_array_finish<Number>(
        vector_operation, temporary_storage, locally_owned_array, requests);
    }



    template <typename Number>
    void
    NoncontiguousPartitioner::import_from_ghosted_array_start(
      const VectorOperation::values vector_operation,
      const unsigned int            communication_channel,
      const ArrayView<Number>      &src,
      const ArrayView<Number>      &buffers,
      std::vector<MPI_Request>     &requests) const
    {
#ifndef DEAL_II_WITH_MPI
      (void)vector_operation;
      (void)communication_channel;
      (void)src;
      (void)buffers;
      (void)requests;
      Assert(false, ExcNeedsMPI());
#else
      (void)vector_operation; // nothing to do here

      AssertDimension(requests.size(), recv_ranks.size() + send_ranks.size());

      const auto tag =
        communication_channel +
        internal::Tags::noncontiguous_partitioner_update_ghost_values_start;

      AssertIndexRange(
        tag,
        internal::Tags::noncontiguous_partitioner_update_ghost_values_end + 1);

      // post recv
      AssertIndexRange(send_ranks.size(), send_ptr.size());
      for (types::global_dof_index i = 0; i < send_ranks.size(); ++i)
        {
          const int ierr =
            MPI_Irecv(buffers.data() + send_ptr[i],
                      send_ptr[i + 1] - send_ptr[i],
                      Utilities::MPI::mpi_type_id_for_type<Number>,
                      send_ranks[i],
                      tag,
                      communicator,
                      &requests[i]);
          AssertThrowMPI(ierr);
        }

      // pack data and send away
      for (types::global_dof_index i = 0; i < recv_ranks.size(); ++i)
        {
          AssertIndexRange(i + 1, recv_ptr.size());
          for (types::global_dof_index j = recv_ptr[i], c = 0;
               j < recv_ptr[i + 1];
               j++)
            buffers[recv_ptr[i] + c++] = src[recv_indices[j]];

          // send data
          Assert((recv_ptr[i] < buffers.size()) ||
                   (recv_ptr[i] == buffers.size() &&
                    recv_ptr[i + 1] == recv_ptr[i]),
                 ExcMessage("The input buffer doesn't contain enough entries"));
          const int ierr =
            MPI_Isend(buffers.data() + recv_ptr[i],
                      recv_ptr[i + 1] - recv_ptr[i],
                      Utilities::MPI::mpi_type_id_for_type<Number>,
                      recv_ranks[i],
                      tag,
                      communicator,
                      &requests[i + send_ranks.size()]);
          AssertThrowMPI(ierr);
        }
#endif
    }



    template <typename Number>
    void
    NoncontiguousPartitioner::import_from_ghosted_array_finish(
      const VectorOperation::values  vector_operation,
      const ArrayView<const Number> &buffers,
      const ArrayView<Number>       &dst,
      std::vector<MPI_Request>      &requests) const
    {
#ifndef DEAL_II_WITH_MPI
      (void)vector_operation;
      (void)buffers;
      (void)dst;
      (void)requests;
      Assert(false, ExcNeedsMPI());
#else
      (void)vector_operation;
      Assert(vector_operation == VectorOperation::add, ExcNotImplemented());

      AssertDimension(requests.size(), recv_ranks.size() + send_ranks.size());

      Assert(std::accumulate(dst.begin(), dst.end(), 0) == 0,
             ExcMessage("The destination vector has to be empty."));

      // wait that all data packages have been received
      // note: this for-loop cold be merged with the next for-loop,
      // however, for this send_indices would be needed to stored
      // rank by rank
      for (types::global_dof_index proc = 0; proc < send_ranks.size(); ++proc)
        {
          int        i;
          MPI_Status status;
          const auto ierr =
            MPI_Waitany(send_ranks.size(), requests.data(), &i, &status);
          AssertThrowMPI(ierr);
        }

      // write data into destination vector
      for (types::global_dof_index i = 0, k = 0; i < send_ranks.size(); ++i)
        {
          // collect data to be send
          for (types::global_dof_index j = send_ptr[i]; j < send_ptr[i + 1];
               j++)
            {
              AssertIndexRange(k, send_indices.size());
              dst[send_indices[k]] += buffers[j];
              ++k;
            }
        }

      // wait that all data packages have been sent
      const int ierr = MPI_Waitall(recv_ranks.size(),
                                   requests.data() + send_ranks.size(),
                                   MPI_STATUSES_IGNORE);
      AssertThrowMPI(ierr);
#endif
    }

  } // namespace MPI
} // namespace Utilities

DEAL_II_NAMESPACE_CLOSE

#endif
