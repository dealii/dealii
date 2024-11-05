// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/mpi.templates.h>
#include <deal.II/base/mpi_consensus_algorithms.h>
#include <deal.II/base/partitioner.h>
#include <deal.II/base/timer.h>

#include <deal.II/matrix_free/vector_data_exchange.h>

#include <boost/serialization/utility.hpp>

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
      PartitionerWrapper::locally_owned_size() const
      {
        return partitioner->locally_owned_size();
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
        const ArrayView<const double>              &locally_owned_array,
        const std::vector<ArrayView<const double>> &shared_arrays,
        const ArrayView<double>                    &ghost_array,
        const ArrayView<double>                    &temporary_storage,
        std::vector<MPI_Request>                   &requests) const
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
        const ArrayView<const double>              &locally_owned_array,
        const std::vector<ArrayView<const double>> &shared_arrays,
        const ArrayView<double>                    &ghost_array,
        std::vector<MPI_Request>                   &requests) const
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
        const ArrayView<const double>              &locally_owned_array,
        const std::vector<ArrayView<const double>> &shared_arrays,
        const ArrayView<double>                    &ghost_array,
        const ArrayView<double>                    &temporary_storage,
        std::vector<MPI_Request>                   &requests) const
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
        const ArrayView<double>                    &locally_owned_storage,
        const std::vector<ArrayView<const double>> &shared_arrays,
        const ArrayView<double>                    &ghost_array,
        const ArrayView<const double>              &temporary_storage,
        std::vector<MPI_Request>                   &requests) const
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
        const ArrayView<const float>              &locally_owned_array,
        const std::vector<ArrayView<const float>> &shared_arrays,
        const ArrayView<float>                    &ghost_array,
        const ArrayView<float>                    &temporary_storage,
        std::vector<MPI_Request>                  &requests) const
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
        const ArrayView<const float>              &locally_owned_array,
        const std::vector<ArrayView<const float>> &shared_arrays,
        const ArrayView<float>                    &ghost_array,
        std::vector<MPI_Request>                  &requests) const
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
        const ArrayView<const float>              &locally_owned_array,
        const std::vector<ArrayView<const float>> &shared_arrays,
        const ArrayView<float>                    &ghost_array,
        const ArrayView<float>                    &temporary_storage,
        std::vector<MPI_Request>                  &requests) const
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
        const ArrayView<float>                    &locally_owned_storage,
        const std::vector<ArrayView<const float>> &shared_arrays,
        const ArrayView<float>                    &ghost_array,
        const ArrayView<const float>              &temporary_storage,
        std::vector<MPI_Request>                  &requests) const
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



      namespace internal
      {
        std::pair<std::vector<unsigned int>,
                  std::vector<std::pair<unsigned int, unsigned int>>>
        compress_to_contiguous_ranges(
          const std::vector<unsigned int> &sm_export_ptr,
          const std::vector<unsigned int> &sm_export_indices)
        {
          std::vector<unsigned int> recv_ptr = {0};
          std::vector<unsigned int> recv_indices;
          std::vector<unsigned int> recv_len;

          for (unsigned int i = 0; i + 1 < sm_export_ptr.size(); ++i)
            {
              if (sm_export_ptr[i] != sm_export_ptr[i + 1])
                {
                  recv_indices.push_back(sm_export_indices[sm_export_ptr[i]]);
                  recv_len.push_back(1);

                  for (unsigned int j = sm_export_ptr[i] + 1;
                       j < sm_export_ptr[i + 1];
                       j++)
                    if (recv_indices.back() + recv_len.back() !=
                        sm_export_indices[j])
                      {
                        recv_indices.push_back(sm_export_indices[j]);
                        recv_len.push_back(1);
                      }
                    else
                      recv_len.back()++;
                }
              recv_ptr.push_back(recv_indices.size());
            }

          std::pair<std::vector<unsigned int>,
                    std::vector<std::pair<unsigned int, unsigned int>>>
            result;

          result.first = recv_ptr;

          for (unsigned int i = 0; i < recv_indices.size(); ++i)
            result.second.emplace_back(recv_indices[i], recv_len[i]);

          return result;
        }

      } // namespace internal



      Full::Full(
        const std::shared_ptr<const Utilities::MPI::Partitioner> &partitioner,
        const MPI_Comm communicator_sm)
        : comm(partitioner->get_mpi_communicator())
        , comm_sm(communicator_sm)
        , n_local_elements(partitioner->locally_owned_range().n_elements())
        , n_ghost_elements(partitioner->ghost_indices().n_elements())
        , n_global_elements(partitioner->locally_owned_range().size())
      {
#ifndef DEAL_II_WITH_MPI
        Assert(false, ExcNeedsMPI());
#else
        if (Utilities::MPI::job_supports_mpi() == false)
          return; // nothing to do in serial case

        const auto &is_locally_owned = partitioner->locally_owned_range();
        const auto &is_locally_ghost = partitioner->ghost_indices();
        const auto &ghost_indices_within_larger_ghost_set =
          partitioner->ghost_indices_within_larger_ghost_set();

        // temporal data structures
        std::vector<unsigned int> n_ghost_indices_in_larger_set_by_remote_rank;

        std::vector<std::array<unsigned int, 3>> ghost_targets_data;

        std::vector<std::array<unsigned int, 3>> import_targets_data;

        std::vector<unsigned int> sm_ghost_ranks;

        std::vector<unsigned int> sm_import_ranks;

        // temporary uncompressed data structures for ghost_indices_subset_data
        std::vector<unsigned int> ghost_indices_subset_data_ptr = {0};
        std::vector<unsigned int> ghost_indices_subset_data_indices;

        // ... for import_indices_data
        std::vector<unsigned int> import_indices_data_ptr = {0};
        std::vector<unsigned int> import_indices_data_indices;

        // ... for sm_export_data
        std::vector<unsigned int> sm_export_data_ptr = {0};
        std::vector<unsigned int> sm_export_data_indices;

        // ... for sm_export_data_this
        std::vector<unsigned int> sm_export_data_this_ptr = {0};
        std::vector<unsigned int> sm_export_data_this_indices;

        // ... for sm_import_data
        std::vector<unsigned int> sm_import_data_ptr = {};
        std::vector<unsigned int> sm_import_data_indices;

        // ... for sm_import_data_this
        std::vector<unsigned int> sm_import_data_this_ptr = {0};
        std::vector<unsigned int> sm_import_data_this_indices;

        // collect ranks of processes of shared-memory domain
        const auto sm_ranks =
          Utilities::MPI::mpi_processes_within_communicator(comm, comm_sm);

        // determine owners of ghost indices and determine requesters
        const auto [owning_ranks_of_ghosts, rank_to_global_indices] =
          Utilities::MPI::compute_index_owner_and_requesters(is_locally_owned,
                                                             is_locally_ghost,
                                                             comm);

        // decompress ghost_indices_within_larger_ghost_set for simpler
        // data access during setup
        std::vector<unsigned int> shifts_indices;
        for (const auto &pair : ghost_indices_within_larger_ghost_set)
          for (unsigned int k = pair.first; k < pair.second; ++k)
            shifts_indices.push_back(k);

        // process ghost indices
        {
          // collect ghost indices according to owning rank
          std::map<unsigned int, std::vector<types::global_dof_index>>
            rank_to_local_indices;

          for (unsigned int i = 0; i < owning_ranks_of_ghosts.size(); ++i)
            rank_to_local_indices[owning_ranks_of_ghosts[i]].push_back(i);

          unsigned int compressed_offset = 0;

          for (const auto &rank_and_local_indices : rank_to_local_indices)
            {
              const auto sm_ranks_ptr = std::find(sm_ranks.begin(),
                                                  sm_ranks.end(),
                                                  rank_and_local_indices.first);

              if (sm_ranks_ptr == sm_ranks.end()) // remote process
                {
                  ghost_targets_data.emplace_back(std::array<unsigned int, 3>{{
                    rank_and_local_indices.first,      // rank
                    shifts_indices[compressed_offset], // offset
                    static_cast<unsigned int>(
                      rank_and_local_indices.second.size()) // length
                  }});

                  for (unsigned int i = 0;
                       i < rank_and_local_indices.second.size();
                       ++i)
                    ghost_indices_subset_data_indices.push_back(
                      shifts_indices[i + compressed_offset]);

                  ghost_indices_subset_data_ptr.push_back(
                    ghost_indices_subset_data_indices.size());

                  ghost_indices_subset_data.first.push_back(compressed_offset);

                  unsigned int i =
                    n_ghost_indices_in_larger_set_by_remote_rank.size();

                  n_ghost_indices_in_larger_set_by_remote_rank.push_back(
                    (shifts_indices[ghost_indices_subset_data.first[i] +
                                    (ghost_targets_data[i][2] - 1)] -
                     shifts_indices[ghost_indices_subset_data.first[i]]) +
                    1);
                }
              else // shared process
                {
                  sm_ghost_ranks.push_back(
                    std::distance(sm_ranks.begin(), sm_ranks_ptr));

                  sm_export_data_ptr.push_back(
                    sm_export_data_ptr.back() +
                    rank_and_local_indices.second.size());

                  for (unsigned int i = compressed_offset;
                       i <
                       rank_and_local_indices.second.size() + compressed_offset;
                       ++i)
                    sm_export_data_this_indices.push_back(
                      shifts_indices[i] + is_locally_owned.n_elements());

                  sm_export_data_this_ptr.push_back(
                    sm_export_data_this_indices.size());
                }
              compressed_offset += rank_and_local_indices.second.size();
            }

          sm_export_data_indices.resize(sm_export_data_ptr.back());
        }

        // process requesters
        {
          for (const auto &rank_and_global_indices : rank_to_global_indices)
            {
              const auto sm_ranks_ptr =
                std::find(sm_ranks.begin(),
                          sm_ranks.end(),
                          rank_and_global_indices.first);

              if (sm_ranks_ptr == sm_ranks.end()) // remote process
                {
                  import_targets_data.emplace_back(std::array<unsigned int, 3>{{
                    rank_and_global_indices.first, // rank
                    static_cast<unsigned int>(
                      import_indices_data_indices.size()), // offset
                    static_cast<unsigned int>(
                      rank_and_global_indices.second.n_elements()) // length
                  }});

                  for (const auto i : rank_and_global_indices.second)
                    import_indices_data_indices.push_back(
                      is_locally_owned.index_within_set(i));

                  import_indices_data_ptr.push_back(
                    import_indices_data_indices.size());
                }
              else // shared process
                {
                  sm_import_ranks.push_back(
                    std::distance(sm_ranks.begin(), sm_ranks_ptr));

                  for (const auto i : rank_and_global_indices.second)
                    sm_import_data_this_indices.push_back(
                      is_locally_owned.index_within_set(i));

                  sm_import_data_this_ptr.push_back(
                    sm_import_data_this_indices.size());
                }
            }

          sm_import_data_ptr = sm_import_data_this_ptr;
          sm_import_data_indices.resize(sm_import_data_this_ptr.back());
        }

        // send sm_export_data_this to sm-neighbor -> sm_import_data
        {
          std::vector<MPI_Request> requests(sm_ghost_ranks.size() +
                                            sm_import_ranks.size());

          for (unsigned int i = 0; i < sm_ghost_ranks.size(); ++i)
            {
              const int ierr = MPI_Isend(sm_export_data_this_indices.data() +
                                           sm_export_data_this_ptr[i],
                                         sm_export_data_this_ptr[i + 1] -
                                           sm_export_data_this_ptr[i],
                                         MPI_UNSIGNED,
                                         sm_ghost_ranks[i],
                                         4,
                                         comm_sm,
                                         requests.data() + i);
              AssertThrowMPI(ierr);
            }

          for (unsigned int i = 0; i < sm_import_ranks.size(); ++i)
            {
              const int ierr =
                MPI_Irecv(sm_import_data_indices.data() + sm_import_data_ptr[i],
                          sm_import_data_ptr[i + 1] - sm_import_data_ptr[i],
                          MPI_UNSIGNED,
                          sm_import_ranks[i],
                          4,
                          comm_sm,
                          requests.data() + sm_ghost_ranks.size() + i);
              AssertThrowMPI(ierr);
            }

          const int ierr =
            MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
          AssertThrowMPI(ierr);
        }

        // send sm_import_data_this to sm-neighbor -> sm_export_data_indices
        {
          std::vector<MPI_Request> requests(sm_import_ranks.size() +
                                            sm_ghost_ranks.size());

          for (unsigned int i = 0; i < sm_import_ranks.size(); ++i)
            {
              const int ierr = MPI_Isend(sm_import_data_this_indices.data() +
                                           sm_import_data_this_ptr[i],
                                         sm_import_data_this_ptr[i + 1] -
                                           sm_import_data_this_ptr[i],
                                         MPI_UNSIGNED,
                                         sm_import_ranks[i],
                                         2,
                                         comm_sm,
                                         requests.data() + i);
              AssertThrowMPI(ierr);
            }

          for (unsigned int i = 0; i < sm_ghost_ranks.size(); ++i)
            {
              const int ierr =
                MPI_Irecv(sm_export_data_indices.data() + sm_export_data_ptr[i],
                          sm_export_data_ptr[i + 1] - sm_export_data_ptr[i],
                          MPI_UNSIGNED,
                          sm_ghost_ranks[i],
                          2,
                          comm_sm,
                          requests.data() + sm_import_ranks.size() + i);
              AssertThrowMPI(ierr);
            }

          const int ierr =
            MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
          AssertThrowMPI(ierr);
        }

        // store data structures and, if needed, compress them
        this->n_ghost_indices_in_larger_set_by_remote_rank =
          n_ghost_indices_in_larger_set_by_remote_rank;

        this->ghost_indices_subset_data =
          internal::compress_to_contiguous_ranges(
            ghost_indices_subset_data_ptr, ghost_indices_subset_data_indices);

        this->ghost_targets_data = ghost_targets_data;

        this->import_targets_data = import_targets_data;

        this->import_indices_data =
          internal::compress_to_contiguous_ranges(import_indices_data_ptr,
                                                  import_indices_data_indices);

        this->sm_ghost_ranks = sm_ghost_ranks;

        this->sm_export_data =
          internal::compress_to_contiguous_ranges(sm_export_data_ptr,
                                                  sm_export_data_indices);

        this->sm_export_data_this =
          internal::compress_to_contiguous_ranges(sm_export_data_this_ptr,
                                                  sm_export_data_this_indices);

        this->sm_import_ranks = sm_import_ranks;

        this->sm_import_data =
          internal::compress_to_contiguous_ranges(sm_import_data_ptr,
                                                  sm_import_data_indices);

        this->sm_import_data_this =
          internal::compress_to_contiguous_ranges(sm_import_data_this_ptr,
                                                  sm_import_data_this_indices);

#endif
      }



      void
      Full::export_to_ghosted_array_start(
        const unsigned int                          communication_channel,
        const ArrayView<const double>              &locally_owned_array,
        const std::vector<ArrayView<const double>> &shared_arrays,
        const ArrayView<double>                    &ghost_array,
        const ArrayView<double>                    &temporary_storage,
        std::vector<MPI_Request>                   &requests) const
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
        const ArrayView<const double>              &locally_owned_array,
        const std::vector<ArrayView<const double>> &shared_arrays,
        const ArrayView<double>                    &ghost_array,
        std::vector<MPI_Request>                   &requests) const
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
        const ArrayView<const double>              &locally_owned_array,
        const std::vector<ArrayView<const double>> &shared_arrays,
        const ArrayView<double>                    &ghost_array,
        const ArrayView<double>                    &temporary_storage,
        std::vector<MPI_Request>                   &requests) const
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
        const ArrayView<double>                    &locally_owned_storage,
        const std::vector<ArrayView<const double>> &shared_arrays,
        const ArrayView<double>                    &ghost_array,
        const ArrayView<const double>              &temporary_storage,
        std::vector<MPI_Request>                   &requests) const
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
        const ArrayView<const float>              &locally_owned_array,
        const std::vector<ArrayView<const float>> &shared_arrays,
        const ArrayView<float>                    &ghost_array,
        const ArrayView<float>                    &temporary_storage,
        std::vector<MPI_Request>                  &requests) const
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
        const ArrayView<const float>              &locally_owned_array,
        const std::vector<ArrayView<const float>> &shared_arrays,
        const ArrayView<float>                    &ghost_array,
        std::vector<MPI_Request>                  &requests) const
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
        const ArrayView<const float>              &locally_owned_array,
        const std::vector<ArrayView<const float>> &shared_arrays,
        const ArrayView<float>                    &ghost_array,
        const ArrayView<float>                    &temporary_storage,
        std::vector<MPI_Request>                  &requests) const
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
        const ArrayView<float>                    &locally_owned_storage,
        const std::vector<ArrayView<const float>> &shared_arrays,
        const ArrayView<float>                    &ghost_array,
        const ArrayView<const float>              &temporary_storage,
        std::vector<MPI_Request>                  &requests) const
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
        const ArrayView<const Number>              &data_this,
        const std::vector<ArrayView<const Number>> &data_others,
        const ArrayView<Number>                    &buffer,
        const ArrayView<Number>                    &temporary_storage,
        std::vector<MPI_Request>                   &requests) const
      {
#ifndef DEAL_II_WITH_MPI
        Assert(false, ExcNeedsMPI());

        (void)communication_channel;
        (void)data_this;
        (void)data_others;
        (void)buffer;
        (void)temporary_storage;
        (void)requests;
#else
        (void)data_others;

        requests.resize(sm_import_ranks.size() + sm_ghost_ranks.size() +
                        ghost_targets_data.size() + import_targets_data.size());

        int dummy;
        // receive a signal that relevant sm neighbors are ready
        for (unsigned int i = 0; i < sm_ghost_ranks.size(); ++i)
          {
            const int ierr =
              MPI_Irecv(&dummy,
                        0,
                        MPI_INT,
                        sm_ghost_ranks[i],
                        communication_channel + 0,
                        comm_sm,
                        requests.data() + sm_import_ranks.size() + i);
            AssertThrowMPI(ierr);
          }

        // signal to all relevant sm neighbors that this process is ready
        for (unsigned int i = 0; i < sm_import_ranks.size(); ++i)
          {
            const int ierr = MPI_Isend(&dummy,
                                       0,
                                       MPI_INT,
                                       sm_import_ranks[i],
                                       communication_channel + 0,
                                       comm_sm,
                                       requests.data() + i);
            AssertThrowMPI(ierr);
          }

        // receive data from remote processes
        for (unsigned int i = 0; i < ghost_targets_data.size(); ++i)
          {
            const unsigned int offset =
              n_ghost_indices_in_larger_set_by_remote_rank[i] -
              ghost_targets_data[i][2];

            const int ierr = MPI_Irecv(
              buffer.data() + ghost_targets_data[i][1] + offset,
              ghost_targets_data[i][2],
              Utilities::MPI::mpi_type_id_for_type<decltype(*buffer.data())>,
              ghost_targets_data[i][0],
              communication_channel + 1,
              comm,
              requests.data() + sm_import_ranks.size() + sm_ghost_ranks.size() +
                i);
            AssertThrowMPI(ierr);
          }

        // send data to remote processes
        for (unsigned int i = 0, k = 0; i < import_targets_data.size(); ++i)
          {
            for (unsigned int j = import_indices_data.first[i];
                 j < import_indices_data.first[i + 1];
                 j++)
              for (unsigned int l = 0; l < import_indices_data.second[j].second;
                   l++, k++)
                temporary_storage[k] =
                  data_this[import_indices_data.second[j].first + l];

            // send data away
            const int ierr = MPI_Isend(
              temporary_storage.data() + import_targets_data[i][1],
              import_targets_data[i][2],
              Utilities::MPI::mpi_type_id_for_type<decltype(*data_this.data())>,
              import_targets_data[i][0],
              communication_channel + 1,
              comm,
              requests.data() + sm_import_ranks.size() + sm_ghost_ranks.size() +
                ghost_targets_data.size() + i);
            AssertThrowMPI(ierr);
          }
#endif
      }



      template <typename Number>
      void
      Full::export_to_ghosted_array_finish_impl(
        const ArrayView<const Number>              &data_this,
        const std::vector<ArrayView<const Number>> &data_others,
        const ArrayView<Number>                    &ghost_array,
        std::vector<MPI_Request>                   &requests) const
      {
        (void)data_this;

#ifndef DEAL_II_WITH_MPI
        Assert(false, ExcNeedsMPI());

        (void)data_others;
        (void)ghost_array;
        (void)requests;
#else

        AssertDimension(requests.size(),
                        sm_import_ranks.size() + sm_ghost_ranks.size() +
                          ghost_targets_data.size() +
                          import_targets_data.size());

        const auto split =
          [&](const unsigned int i) -> std::pair<unsigned int, unsigned int> {
          AssertIndexRange(i,
                           (sm_ghost_ranks.size() + ghost_targets_data.size()));

          if (i < sm_ghost_ranks.size())
            return {0, i};
          else
            return {1, i - sm_ghost_ranks.size()};
        };

        for (unsigned int c = 0;
             c < sm_ghost_ranks.size() + ghost_targets_data.size();
             c++)
          {
            int       i;
            const int ierr =
              MPI_Waitany(sm_ghost_ranks.size() + ghost_targets_data.size(),
                          requests.data() + sm_import_ranks.size(),
                          &i,
                          MPI_STATUS_IGNORE);
            AssertThrowMPI(ierr);

            const auto s = split(i);
            i            = s.second;

            if (s.first == 0)
              {
                const Number *DEAL_II_RESTRICT data_others_ptr =
                  data_others[sm_ghost_ranks[i]].data();
                Number *DEAL_II_RESTRICT data_this_ptr = ghost_array.data();

                for (unsigned int lo = sm_export_data.first[i],
                                  ko = sm_export_data_this.first[i],
                                  li = 0,
                                  ki = 0;
                     (lo < sm_export_data.first[i + 1]) &&
                     (ko < sm_export_data_this.first[i + 1]);)
                  {
                    for (; (li < sm_export_data.second[lo].second) &&
                           (ki < sm_export_data_this.second[ko].second);
                         ++li, ++ki)
                      data_this_ptr[sm_export_data_this.second[ko].first + ki -
                                    n_local_elements] =
                        data_others_ptr[sm_export_data.second[lo].first + li];

                    if (li == sm_export_data.second[lo].second)
                      {
                        lo++;   // increment outer counter
                        li = 0; // reset inner counter
                      }

                    if (ki == sm_export_data_this.second[ko].second)
                      {
                        ko++;   // increment outer counter
                        ki = 0; // reset inner counter
                      }
                  }
              }
            else /*if(s.second == 1)*/
              {
                const unsigned int offset =
                  n_ghost_indices_in_larger_set_by_remote_rank[i] -
                  ghost_targets_data[i][2];

                for (unsigned int c  = 0,
                                  ko = ghost_indices_subset_data.first[i],
                                  ki = 0;
                     c < ghost_targets_data[i][2];
                     ++c)
                  {
                    AssertIndexRange(ko,
                                     ghost_indices_subset_data.second.size());

                    const unsigned int idx_1 =
                      ghost_indices_subset_data.second[ko].first + ki;
                    const unsigned int idx_2 =
                      ghost_targets_data[i][1] + c + offset;

                    AssertIndexRange(idx_1, ghost_array.size());
                    AssertIndexRange(idx_2, ghost_array.size());

                    if (idx_1 == idx_2)
                      {
                        // noting to do
                      }
                    else if (idx_1 < idx_2)
                      {
                        ghost_array[idx_1] = ghost_array[idx_2];
                        ghost_array[idx_2] = 0.0;
                      }
                    else
                      {
                        DEAL_II_NOT_IMPLEMENTED();
                      }

                    ++ki;

                    if (ki == ghost_indices_subset_data.second[ko].second)
                      {
                        ko++;   // increment outer counter
                        ki = 0; // reset inner counter
                      }
                  }
              }
          }

        const int ierr =
          MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
        AssertThrowMPI(ierr);

#endif
      }



      template <typename Number>
      void
      Full::import_from_ghosted_array_start_impl(
        const VectorOperation::values               operation,
        const unsigned int                          communication_channel,
        const ArrayView<const Number>              &data_this,
        const std::vector<ArrayView<const Number>> &data_others,
        const ArrayView<Number>                    &buffer,
        const ArrayView<Number>                    &temporary_storage,
        std::vector<MPI_Request>                   &requests) const
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
        // return;

        (void)data_others;
        (void)operation;

        Assert(operation == VectorOperation::add, ExcNotImplemented());

        requests.resize(sm_ghost_ranks.size() + sm_import_ranks.size() +
                        ghost_targets_data.size() + import_targets_data.size());

        int dummy;
        for (unsigned int i = 0; i < sm_ghost_ranks.size(); ++i)
          {
            const int ierr = MPI_Isend(&dummy,
                                       0,
                                       MPI_INT,
                                       sm_ghost_ranks[i],
                                       communication_channel + 1,
                                       comm_sm,
                                       requests.data() + i);
            AssertThrowMPI(ierr);
          }

        for (unsigned int i = 0; i < sm_import_ranks.size(); ++i)
          {
            const int ierr =
              MPI_Irecv(&dummy,
                        0,
                        MPI_INT,
                        sm_import_ranks[i],
                        communication_channel + 1,
                        comm_sm,
                        requests.data() + sm_ghost_ranks.size() + i);
            AssertThrowMPI(ierr);
          }

        for (unsigned int i = 0; i < ghost_targets_data.size(); ++i)
          {
            for (unsigned int c  = 0,
                              ko = ghost_indices_subset_data.first[i],
                              ki = 0;
                 c < ghost_targets_data[i][2];
                 ++c)
              {
                AssertIndexRange(ko, ghost_indices_subset_data.second.size());

                const unsigned int idx_1 =
                  ghost_indices_subset_data.second[ko].first + ki;
                const unsigned int idx_2 = ghost_targets_data[i][1] + c;

                AssertIndexRange(idx_1, buffer.size());
                AssertIndexRange(idx_2, buffer.size());

                if (idx_1 == idx_2)
                  {
                    // nothing to do
                  }
                else if (idx_2 < idx_1)
                  {
                    buffer[idx_2] = buffer[idx_1];
                    buffer[idx_1] = 0.0;
                  }
                else
                  {
                    DEAL_II_NOT_IMPLEMENTED();
                  }

                if (++ki == ghost_indices_subset_data.second[ko].second)
                  {
                    ko++;   // increment outer counter
                    ki = 0; // reset inner counter
                  }
              }

            const int ierr = MPI_Isend(
              buffer.data() + ghost_targets_data[i][1],
              ghost_targets_data[i][2],
              Utilities::MPI::mpi_type_id_for_type<decltype(*buffer.data())>,
              ghost_targets_data[i][0],
              communication_channel + 0,
              comm,
              requests.data() + sm_ghost_ranks.size() + sm_import_ranks.size() +
                i);
            AssertThrowMPI(ierr);
          }

        for (unsigned int i = 0; i < import_targets_data.size(); ++i)
          {
            const int ierr =
              MPI_Irecv(temporary_storage.data() + import_targets_data[i][1],
                        import_targets_data[i][2],
                        Utilities::MPI::mpi_type_id_for_type<
                          decltype(*temporary_storage.data())>,
                        import_targets_data[i][0],
                        communication_channel + 0,
                        comm,
                        requests.data() + sm_ghost_ranks.size() +
                          sm_import_ranks.size() + ghost_targets_data.size() +
                          i);
            AssertThrowMPI(ierr);
          }
#endif
      }



      template <typename Number>
      void
      Full::import_from_ghosted_array_finish_impl(
        const VectorOperation::values               operation,
        const ArrayView<Number>                    &data_this,
        const std::vector<ArrayView<const Number>> &data_others,
        const ArrayView<Number>                    &buffer,
        const ArrayView<const Number>              &temporary_storage,
        std::vector<MPI_Request>                   &requests) const
      {
#ifndef DEAL_II_WITH_MPI
        Assert(false, ExcNeedsMPI());

        (void)operation;
        (void)data_this;
        (void)data_others;
        (void)buffer;
        (void)temporary_storage;
        (void)requests;
#else

        (void)operation;

        Assert(operation == VectorOperation::add, ExcNotImplemented());

        AssertDimension(requests.size(),
                        sm_ghost_ranks.size() + sm_import_ranks.size() +
                          ghost_targets_data.size() +
                          import_targets_data.size());

        const auto split =
          [&](const unsigned int i) -> std::pair<unsigned int, unsigned int> {
          AssertIndexRange(i,
                           (sm_import_ranks.size() + ghost_targets_data.size() +
                            import_targets_data.size()));

          if (i < sm_import_ranks.size())
            return {0, i};
          else if (i < (sm_import_ranks.size() + ghost_targets_data.size()))
            return {2, i - sm_import_ranks.size()};
          else
            return {1, i - sm_import_ranks.size() - ghost_targets_data.size()};
        };

        for (unsigned int c = 0;
             c < sm_import_ranks.size() + import_targets_data.size() +
                   ghost_targets_data.size();
             c++)
          {
            int       i;
            const int ierr =
              MPI_Waitany(sm_import_ranks.size() + import_targets_data.size() +
                            ghost_targets_data.size(),
                          requests.data() + sm_ghost_ranks.size(),
                          &i,
                          MPI_STATUS_IGNORE);
            AssertThrowMPI(ierr);

            const auto &s = split(i);
            i             = s.second;

            if (s.first == 0)
              {
                Number *DEAL_II_RESTRICT data_others_ptr =
                  const_cast<Number *>(data_others[sm_import_ranks[i]].data());
                Number *DEAL_II_RESTRICT data_this_ptr = data_this.data();

                for (unsigned int lo = sm_import_data_this.first[i],
                                  ko = sm_import_data.first[i],
                                  li = 0,
                                  ki = 0;
                     (lo < sm_import_data_this.first[i + 1]) &&
                     (ko < sm_import_data.first[i + 1]);)
                  {
                    for (; (li < sm_import_data_this.second[lo].second) &&
                           (ki < sm_import_data.second[ko].second);
                         ++li, ++ki)
                      {
                        data_this_ptr[sm_import_data_this.second[lo].first +
                                      li] +=
                          data_others_ptr[sm_import_data.second[ko].first + ki];
                        data_others_ptr[sm_import_data.second[ko].first + ki] =
                          0.0;
                      }

                    if (li == sm_import_data_this.second[lo].second)
                      {
                        lo++;   // increment outer counter
                        li = 0; // reset inner counter
                      }
                    if (ki == sm_import_data.second[ko].second)
                      {
                        ko++;   // increment outer counter
                        ki = 0; // reset inner counter
                      }
                  }
              }
            else if (s.first == 1)
              {
                for (unsigned int j = import_indices_data.first[i],
                                  k = import_targets_data[i][1];
                     j < import_indices_data.first[i + 1];
                     j++)
                  for (unsigned int l = 0;
                       l < import_indices_data.second[j].second;
                       l++)
                    data_this[import_indices_data.second[j].first + l] +=
                      temporary_storage[k++];
              }
            else /*if (s.first == 2)*/
              {
                std::memset(buffer.data() + ghost_targets_data[i][1],
                            0.0,
                            (ghost_targets_data[i][2]) * sizeof(Number));
              }
          }

        const int ierr =
          MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
        AssertThrowMPI(ierr);
#endif
      }



      unsigned int
      Full::locally_owned_size() const
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
        if (import_targets_data.empty())
          return 0;
        return import_targets_data.back()[1] + import_targets_data.back()[2];
      }



      unsigned int
      Full::n_import_sm_procs() const
      {
        return sm_import_ranks.size() + sm_ghost_ranks.size(); // TODO
      }



      types::global_dof_index
      Full::size() const
      {
        return n_global_elements;
      }



      MPI_Comm
      Full::get_sm_mpi_communicator() const
      {
        return this->comm_sm;
      }



      void
      Full::reset_ghost_values(const ArrayView<double> &ghost_array) const
      {
        reset_ghost_values_impl(ghost_array);
      }



      void
      Full::reset_ghost_values(const ArrayView<float> &ghost_array) const
      {
        reset_ghost_values_impl(ghost_array);
      }



      template <typename Number>
      void
      Full::reset_ghost_values_impl(const ArrayView<Number> &ghost_array) const
      {
        // reset ghost values coming from shared-memory neighbors
        // TODO: only needed if values are buffered
        for (const auto &i : sm_export_data_this.second)
          std::memset(ghost_array.data() + (i.first - n_local_elements),
                      0,
                      sizeof(Number) * i.second);

        // reset ghost values coming from remote neighbors
        for (const auto &i : ghost_indices_subset_data.second)
          std::memset(ghost_array.data() + i.first,
                      0,
                      sizeof(Number) * i.second);
      }



    } // namespace VectorDataExchange
  }   // namespace MatrixFreeFunctions
} // namespace internal


DEAL_II_NAMESPACE_CLOSE
