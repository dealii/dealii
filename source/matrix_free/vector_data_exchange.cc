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

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/mpi.templates.h>
#include <deal.II/base/mpi_compute_index_owner_internal.h>
#include <deal.II/base/mpi_consensus_algorithms.h>
#include <deal.II/base/timer.h>

#include <deal.II/matrix_free/vector_data_exchange.h>

#ifdef DEAL_II_WITH_64BIT_INDICES
#  include <deal.II/base/mpi_consensus_algorithms.templates.h>
#endif

#include <map>
#include <vector>


DEAL_II_NAMESPACE_OPEN

namespace internal
{
  namespace MatrixFreeFunctions
  {
    namespace VectorDataExchange
    {
      PartitionerWrapper::PartitionerWrapper(
        const std::shared_ptr<const Utilities::MPI::Partitioner> &partitioner)
        : partitioner(partitioner)
      {}



      unsigned int
      PartitionerWrapper::local_size() const
      {
        return partitioner->local_size();
      }



      unsigned int
      PartitionerWrapper::n_ghost_indices() const
      {
        return partitioner->n_ghost_indices();
      }



      unsigned int
      PartitionerWrapper::n_import_indices() const
      {
        return partitioner->n_import_indices();
      }



      unsigned int
      PartitionerWrapper::n_import_sm_procs() const
      {
        return 0;
      }



      types::global_dof_index
      PartitionerWrapper::size() const
      {
        return partitioner->size();
      }



      void
      PartitionerWrapper::export_to_ghosted_array_start(
        const unsigned int                          communication_channel,
        const ArrayView<const double> &             locally_owned_array,
        const std::vector<ArrayView<const double>> &shared_arrays,
        const ArrayView<double> &                   ghost_array,
        const ArrayView<double> &                   temporary_storage,
        std::vector<MPI_Request> &                  requests) const
      {
        (void)shared_arrays;
#ifndef DEAL_II_WITH_MPI
        (void)communication_channel;
        (void)locally_owned_array;
        (void)ghost_array;
        (void)temporary_storage;
        (void)requests;
#else
        partitioner->export_to_ghosted_array_start(communication_channel,
                                                   locally_owned_array,
                                                   temporary_storage,
                                                   ghost_array,
                                                   requests);
#endif
      }



      void
      PartitionerWrapper::export_to_ghosted_array_finish(
        const ArrayView<const double> &             locally_owned_array,
        const std::vector<ArrayView<const double>> &shared_arrays,
        const ArrayView<double> &                   ghost_array,
        std::vector<MPI_Request> &                  requests) const
      {
        (void)locally_owned_array;
        (void)shared_arrays;
#ifndef DEAL_II_WITH_MPI
        (void)ghost_array;
        (void)requests;
#else
        partitioner->export_to_ghosted_array_finish(ghost_array, requests);
#endif
      }



      void
      PartitionerWrapper::import_from_ghosted_array_start(
        const VectorOperation::values               vector_operation,
        const unsigned int                          communication_channel,
        const ArrayView<const double> &             locally_owned_array,
        const std::vector<ArrayView<const double>> &shared_arrays,
        const ArrayView<double> &                   ghost_array,
        const ArrayView<double> &                   temporary_storage,
        std::vector<MPI_Request> &                  requests) const
      {
        (void)locally_owned_array;
        (void)shared_arrays;
#ifndef DEAL_II_WITH_MPI
        (void)vector_operation;
        (void)communication_channel;
        (void)ghost_array;
        (void)temporary_storage;
        (void)requests;
#else
        partitioner->import_from_ghosted_array_start(vector_operation,
                                                     communication_channel,
                                                     ghost_array,
                                                     temporary_storage,
                                                     requests);
#endif
      }



      void
      PartitionerWrapper::import_from_ghosted_array_finish(
        const VectorOperation::values               vector_operation,
        const ArrayView<double> &                   locally_owned_storage,
        const std::vector<ArrayView<const double>> &shared_arrays,
        const ArrayView<double> &                   ghost_array,
        const ArrayView<const double> &             temporary_storage,
        std::vector<MPI_Request> &                  requests) const
      {
        (void)shared_arrays;
#ifndef DEAL_II_WITH_MPI
        (void)vector_operation;
        (void)locally_owned_storage;
        (void)ghost_array;
        (void)temporary_storage;
        (void)requests;
#else
        partitioner->import_from_ghosted_array_finish(vector_operation,
                                                      temporary_storage,
                                                      locally_owned_storage,
                                                      ghost_array,
                                                      requests);
#endif
      }



      void
      PartitionerWrapper::reset_ghost_values(
        const ArrayView<double> &ghost_array) const
      {
        reset_ghost_values_impl(ghost_array);
      }



      void
      PartitionerWrapper::export_to_ghosted_array_start(
        const unsigned int                         communication_channel,
        const ArrayView<const float> &             locally_owned_array,
        const std::vector<ArrayView<const float>> &shared_arrays,
        const ArrayView<float> &                   ghost_array,
        const ArrayView<float> &                   temporary_storage,
        std::vector<MPI_Request> &                 requests) const
      {
        (void)shared_arrays;
#ifndef DEAL_II_WITH_MPI
        (void)communication_channel;
        (void)locally_owned_array;
        (void)ghost_array;
        (void)temporary_storage;
        (void)requests;
#else
        partitioner->export_to_ghosted_array_start(communication_channel,
                                                   locally_owned_array,
                                                   temporary_storage,
                                                   ghost_array,
                                                   requests);
#endif
      }



      void
      PartitionerWrapper::export_to_ghosted_array_finish(
        const ArrayView<const float> &             locally_owned_array,
        const std::vector<ArrayView<const float>> &shared_arrays,
        const ArrayView<float> &                   ghost_array,
        std::vector<MPI_Request> &                 requests) const
      {
        (void)locally_owned_array;
        (void)shared_arrays;
#ifndef DEAL_II_WITH_MPI
        (void)ghost_array;
        (void)requests;
#else
        partitioner->export_to_ghosted_array_finish(ghost_array, requests);
#endif
      }



      void
      PartitionerWrapper::import_from_ghosted_array_start(
        const VectorOperation::values              vector_operation,
        const unsigned int                         communication_channel,
        const ArrayView<const float> &             locally_owned_array,
        const std::vector<ArrayView<const float>> &shared_arrays,
        const ArrayView<float> &                   ghost_array,
        const ArrayView<float> &                   temporary_storage,
        std::vector<MPI_Request> &                 requests) const
      {
        (void)locally_owned_array;
        (void)shared_arrays;
#ifndef DEAL_II_WITH_MPI
        (void)vector_operation;
        (void)communication_channel;
        (void)ghost_array;
        (void)temporary_storage;
        (void)requests;
#else
        partitioner->import_from_ghosted_array_start(vector_operation,
                                                     communication_channel,
                                                     ghost_array,
                                                     temporary_storage,
                                                     requests);
#endif
      }



      void
      PartitionerWrapper::import_from_ghosted_array_finish(
        const VectorOperation::values              vector_operation,
        const ArrayView<float> &                   locally_owned_storage,
        const std::vector<ArrayView<const float>> &shared_arrays,
        const ArrayView<float> &                   ghost_array,
        const ArrayView<const float> &             temporary_storage,
        std::vector<MPI_Request> &                 requests) const
      {
        (void)shared_arrays;
#ifndef DEAL_II_WITH_MPI
        (void)vector_operation;
        (void)locally_owned_storage;
        (void)ghost_array;
        (void)temporary_storage;
        (void)requests;
#else
        partitioner->import_from_ghosted_array_finish(vector_operation,
                                                      temporary_storage,
                                                      locally_owned_storage,
                                                      ghost_array,
                                                      requests);
#endif
      }



      void
      PartitionerWrapper::reset_ghost_values(
        const ArrayView<float> &ghost_array) const
      {
        reset_ghost_values_impl(ghost_array);
      }

      template <typename Number>
      void
      PartitionerWrapper::reset_ghost_values_impl(
        const ArrayView<Number> &ghost_array) const
      {
        for (const auto &my_ghosts :
             partitioner->ghost_indices_within_larger_ghost_set())
          for (unsigned int j = my_ghosts.first; j < my_ghosts.second; ++j)
            ghost_array[j] = 0.;
      }

    } // namespace VectorDataExchange
  }   // namespace MatrixFreeFunctions
} // namespace internal


DEAL_II_NAMESPACE_CLOSE
