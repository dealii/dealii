// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_partitioner_templates_h
#define dealii_partitioner_templates_h

#include <deal.II/base/config.h>

#include <deal.II/base/mpi_tags.h>
#include <deal.II/base/partitioner.h>

#include <deal.II/lac/la_parallel_vector.h>

#include <limits>
#include <type_traits>


DEAL_II_NAMESPACE_OPEN

namespace Utilities
{
  namespace MPI
  {
#ifndef DOXYGEN

#  ifdef DEAL_II_WITH_MPI

    template <typename Number, typename MemorySpaceType>
    void
    Partitioner::export_to_ghosted_array_start(
      const unsigned int                              communication_channel,
      const ArrayView<const Number, MemorySpaceType> &locally_owned_array,
      const ArrayView<Number, MemorySpaceType>       &temporary_storage,
      const ArrayView<Number, MemorySpaceType>       &ghost_array,
      std::vector<MPI_Request>                       &requests) const
    {
      AssertDimension(temporary_storage.size(), n_import_indices());
      AssertIndexRange(communication_channel, 200);
      Assert(ghost_array.size() == n_ghost_indices() ||
               ghost_array.size() == n_ghost_indices_in_larger_set,
             ExcGhostIndexArrayHasWrongSize(ghost_array.size(),
                                            n_ghost_indices(),
                                            n_ghost_indices_in_larger_set));

      const unsigned int n_import_targets = import_targets_data.size();
      const unsigned int n_ghost_targets  = ghost_targets_data.size();

      if (n_import_targets > 0)
        AssertDimension(locally_owned_array.size(), locally_owned_size());

      Assert(requests.empty(),
             ExcMessage("Another operation seems to still be running. "
                        "Call update_ghost_values_finish() first."));

      const unsigned int mpi_tag =
        Utilities::MPI::internal::Tags::partitioner_export_start +
        communication_channel;
      Assert(mpi_tag <= Utilities::MPI::internal::Tags::partitioner_export_end,
             ExcInternalError());

      // Need to send and receive the data. Use non-blocking communication,
      // where it is usually less overhead to first initiate the receive and
      // then actually send the data
      requests.resize(n_import_targets + n_ghost_targets);

      // as a ghost array pointer, put the data at the end of the given ghost
      // array in case we want to fill only a subset of the ghosts so that we
      // can move data to the right position in a forward loop in the _finish
      // function.
      AssertIndexRange(n_ghost_indices(), n_ghost_indices_in_larger_set + 1);
      const bool use_larger_set =
        (n_ghost_indices_in_larger_set > n_ghost_indices() &&
         ghost_array.size() == n_ghost_indices_in_larger_set);
      Number *ghost_array_ptr =
        use_larger_set ? ghost_array.data() + n_ghost_indices_in_larger_set -
                           n_ghost_indices() :
                         ghost_array.data();

      for (unsigned int i = 0; i < n_ghost_targets; ++i)
        {
          // allow writing into ghost indices even though we are in a
          // const function
          const int ierr =
            MPI_Irecv(ghost_array_ptr,
                      ghost_targets_data[i].second * sizeof(Number),
                      MPI_BYTE,
                      ghost_targets_data[i].first,
                      mpi_tag,
                      communicator,
                      &requests[i]);
          AssertThrowMPI(ierr);
          ghost_array_ptr += ghost_targets_data[i].second;
        }

      Number *temp_array_ptr = temporary_storage.data();
#    if defined(DEAL_II_MPI_WITH_DEVICE_SUPPORT)
      // When using device-aware MPI, the set of local indices that are ghosts
      // indices on other processors is expanded in arrays. This is for
      // performance reasons as this can significantly decrease the number of
      // kernel launched. The indices are expanded the first time the function
      // is called.
      if ((std::is_same_v<MemorySpaceType, MemorySpace::Default>)&&(
            import_indices_plain_dev.empty()))
        initialize_import_indices_plain_dev();
#    endif

      for (unsigned int i = 0; i < n_import_targets; ++i)
        {
#    if defined(DEAL_II_MPI_WITH_DEVICE_SUPPORT)
          if constexpr (std::is_same_v<MemorySpaceType, MemorySpace::Default>)
            {
              const auto chunk_size = import_indices_plain_dev[i].size();
              using IndexType       = decltype(chunk_size);

              auto import_indices           = import_indices_plain_dev[i];
              auto locally_owned_array_data = locally_owned_array.data();
              MemorySpace::Default::kokkos_space::execution_space exec;
              Kokkos::parallel_for(
                "dealii::fill temp_array_ptr",
                Kokkos::RangePolicy<
                  MemorySpace::Default::kokkos_space::execution_space>(
                  exec, 0, chunk_size),
                KOKKOS_LAMBDA(IndexType idx) {
                  temp_array_ptr[idx] =
                    locally_owned_array_data[import_indices[idx]];
                });
              exec.fence();
            }
          else
#    endif
            {
              // copy the data to be sent to the import_data field
              std::vector<std::pair<unsigned int, unsigned int>>::const_iterator
                my_imports = import_indices_data.begin() +
                             import_indices_chunks_by_rank_data[i],
                end_my_imports = import_indices_data.begin() +
                                 import_indices_chunks_by_rank_data[i + 1];
              unsigned int index = 0;
              for (; my_imports != end_my_imports; ++my_imports)
                {
                  const unsigned int chunk_size =
                    my_imports->second - my_imports->first;
                  {
                    std::memcpy(temp_array_ptr + index,
                                locally_owned_array.data() + my_imports->first,
                                chunk_size * sizeof(Number));
                  }
                  index += chunk_size;
                }

              AssertDimension(index, import_targets_data[i].second);
            }

          // start the send operations
          const int ierr =
            MPI_Isend(temp_array_ptr,
                      import_targets_data[i].second * sizeof(Number),
                      MPI_BYTE,
                      import_targets_data[i].first,
                      mpi_tag,
                      communicator,
                      &requests[n_ghost_targets + i]);
          AssertThrowMPI(ierr);
          temp_array_ptr += import_targets_data[i].second;
        }
    }



    template <typename Number, typename MemorySpaceType>
    void
    Partitioner::export_to_ghosted_array_finish(
      const ArrayView<Number, MemorySpaceType> &ghost_array,
      std::vector<MPI_Request>                 &requests) const
    {
      Assert(ghost_array.size() == n_ghost_indices() ||
               ghost_array.size() == n_ghost_indices_in_larger_set,
             ExcGhostIndexArrayHasWrongSize(ghost_array.size(),
                                            n_ghost_indices(),
                                            n_ghost_indices_in_larger_set));

      // wait for both sends and receives to complete, even though only
      // receives are really necessary. this gives (much) better performance
      AssertDimension(ghost_targets().size() + import_targets().size(),
                      requests.size());
      if (requests.size() > 0)
        {
          const int ierr =
            MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
          AssertThrowMPI(ierr);
        }
      requests.resize(0);

      // in case we only sent a subset of indices, we now need to move the data
      // to the correct positions and delete the old content
      if (n_ghost_indices_in_larger_set > n_ghost_indices() &&
          ghost_array.size() == n_ghost_indices_in_larger_set)
        {
          unsigned int offset =
            n_ghost_indices_in_larger_set - n_ghost_indices();
          // must copy ghost data into extended ghost array
          for (const auto &ghost_range : ghost_indices_subset_data)
            {
              if (offset > ghost_range.first)
                {
                  const unsigned int chunk_size =
                    ghost_range.second - ghost_range.first;
                  if constexpr (std::is_same_v<MemorySpaceType,
                                               MemorySpace::Host>)
                    {
                      // If source and destination are overlapping, we must be
                      // careful to use an appropriate copy function.
                      if (ghost_range.first > offset)
                        std::copy_backward(ghost_array.data() + offset,
                                           ghost_array.data() + offset +
                                             chunk_size,
                                           ghost_array.data() +
                                             ghost_range.first);
                      else
                        std::copy(ghost_array.data() + offset,
                                  ghost_array.data() + offset + chunk_size,
                                  ghost_array.data() + ghost_range.first);
                      std::fill(ghost_array.data() +
                                  std::max(ghost_range.second, offset),
                                ghost_array.data() + offset + chunk_size,
                                Number{});
                    }
                  else
                    {
                      Kokkos::View<const Number *,
                                   MemorySpace::Default::kokkos_space>
                        ghost_src_view(ghost_array.data() + offset, chunk_size);
                      Kokkos::View<Number *, MemorySpace::Default::kokkos_space>
                        ghost_dst_view(ghost_array.data() + ghost_range.first,
                                       chunk_size);

                      // If source and destination are overlapping, we can't
                      // just call deep_copy but must copy the data to a buffer
                      // first.
                      if ((offset < ghost_range.first &&
                           ghost_range.first < offset + chunk_size) ||
                          (ghost_range.first < offset &&
                           offset < ghost_range.first + chunk_size))
                        {
                          Kokkos::View<Number *,
                                       MemorySpace::Default::kokkos_space>
                            copy(Kokkos::view_alloc(
                                   "copy", Kokkos::WithoutInitializing),
                                 chunk_size);
                          Kokkos::deep_copy(copy, ghost_src_view);
                          Kokkos::deep_copy(ghost_dst_view, copy);
                        }
                      else
                        {
                          Kokkos::deep_copy(ghost_dst_view, ghost_src_view);
                        }
                      Kokkos::deep_copy(
                        Kokkos::View<Number *,
                                     MemorySpace::Default::kokkos_space>(
                          ghost_array.data() +
                            std::max(ghost_range.second, offset),
                          (offset + chunk_size -
                           std::max(ghost_range.second, offset))),
                        0);
                    }
                  offset += chunk_size;
                }
              else
                {
                  AssertDimension(offset, ghost_range.first);
                  break;
                }
            }
        }
    }



    template <typename Number, typename MemorySpaceType>
    void
    Partitioner::import_from_ghosted_array_start(
      const VectorOperation::values             vector_operation,
      const unsigned int                        communication_channel,
      const ArrayView<Number, MemorySpaceType> &ghost_array,
      const ArrayView<Number, MemorySpaceType> &temporary_storage,
      std::vector<MPI_Request>                 &requests) const
    {
      AssertDimension(temporary_storage.size(), n_import_indices());
      AssertIndexRange(communication_channel, 200);
      Assert(ghost_array.size() == n_ghost_indices() ||
               ghost_array.size() == n_ghost_indices_in_larger_set,
             ExcGhostIndexArrayHasWrongSize(ghost_array.size(),
                                            n_ghost_indices(),
                                            n_ghost_indices_in_larger_set));

      (void)vector_operation;

      // nothing to do for insert (only need to zero ghost entries in
      // compress_finish()). in debug mode we want to check consistency of the
      // inserted data, therefore the communication is still initialized.
      // Having different code in debug and optimized mode is somewhat
      // dangerous, but it really saves communication so do it anyway
      if constexpr (running_in_debug_mode() == false)
        if (vector_operation == VectorOperation::insert)
          return;

      // nothing to do when we neither have import
      // nor ghost indices.
      if (n_ghost_indices() == 0 && n_import_indices() == 0)
        return;

      const unsigned int n_import_targets = import_targets_data.size();
      const unsigned int n_ghost_targets  = ghost_targets_data.size();

      Assert(requests.empty(),
             ExcMessage("Another compress operation seems to still be running. "
                        "Call compress_finish() first."));

      // Need to send and receive the data. Use non-blocking communication,
      // where it is generally less overhead to first initiate the receive and
      // then actually send the data

      const unsigned int mpi_tag =
        Utilities::MPI::internal::Tags::partitioner_import_start +
        communication_channel;
      Assert(mpi_tag <= Utilities::MPI::internal::Tags::partitioner_import_end,
             ExcInternalError());
      requests.resize(n_import_targets + n_ghost_targets);

      // initiate the receive operations
      Number *temp_array_ptr = temporary_storage.data();
      for (unsigned int i = 0; i < n_import_targets; ++i)
        {
          AssertThrow(
            static_cast<std::size_t>(import_targets_data[i].second) *
                sizeof(Number) <
              static_cast<std::size_t>(std::numeric_limits<int>::max()),
            ExcMessage("Index overflow: Maximum message size in MPI is 2GB. "
                       "The number of ghost entries times the size of 'Number' "
                       "exceeds this value. This is not supported."));
          const int ierr =
            MPI_Irecv(temp_array_ptr,
                      import_targets_data[i].second * sizeof(Number),
                      MPI_BYTE,
                      import_targets_data[i].first,
                      mpi_tag,
                      communicator,
                      &requests[i]);
          AssertThrowMPI(ierr);
          temp_array_ptr += import_targets_data[i].second;
        }

      // initiate the send operations

      // in case we want to import only from a subset of the ghosts we want to
      // move the data to send to the front of the array
      AssertIndexRange(n_ghost_indices(), n_ghost_indices_in_larger_set + 1);
      Number *ghost_array_ptr = ghost_array.data();
      for (unsigned int i = 0; i < n_ghost_targets; ++i)
        {
          // in case we only sent a subset of indices, we now need to move the
          // data to the correct positions and delete the old content
          if (n_ghost_indices_in_larger_set > n_ghost_indices() &&
              ghost_array.size() == n_ghost_indices_in_larger_set)
            {
              std::vector<std::pair<unsigned int, unsigned int>>::const_iterator
                my_ghosts = ghost_indices_subset_data.begin() +
                            ghost_indices_subset_chunks_by_rank_data[i],
                end_my_ghosts = ghost_indices_subset_data.begin() +
                                ghost_indices_subset_chunks_by_rank_data[i + 1];
              unsigned int offset = 0;
              for (; my_ghosts != end_my_ghosts; ++my_ghosts)
                {
                  const unsigned int chunk_size =
                    my_ghosts->second - my_ghosts->first;
                  if (ghost_array_ptr + offset !=
                      ghost_array.data() + my_ghosts->first)
                    {
                      if constexpr (std::is_same_v<MemorySpaceType,
                                                   MemorySpace::Host>)
                        {
                          if (offset > my_ghosts->first)
                            std::copy_backward(ghost_array.data() +
                                                 my_ghosts->first,
                                               ghost_array_ptr +
                                                 my_ghosts->second,
                                               ghost_array.data() + offset);
                          else
                            std::copy(ghost_array.data() + my_ghosts->first,
                                      ghost_array.data() + my_ghosts->second,
                                      ghost_array_ptr + offset);
                          std::fill(
                            std::max(ghost_array.data() + my_ghosts->first,
                                     ghost_array_ptr + offset + chunk_size),
                            ghost_array.data() + my_ghosts->second,
                            Number{});
                        }
                      else
                        {
                          Kokkos::View<Number *,
                                       MemorySpace::Default::kokkos_space>
                            copy("copy", chunk_size);
                          Kokkos::deep_copy(
                            copy,
                            Kokkos::View<Number *,
                                         MemorySpace::Default::kokkos_space>(
                              ghost_array.data() + my_ghosts->first,
                              chunk_size));
                          Kokkos::deep_copy(
                            Kokkos::View<Number *,
                                         MemorySpace::Default::kokkos_space>(
                              ghost_array_ptr + offset, chunk_size),
                            copy);
                          Kokkos::deep_copy(
                            Kokkos::View<Number *,
                                         MemorySpace::Default::kokkos_space>(
                              std::max(ghost_array.data() + my_ghosts->first,
                                       ghost_array_ptr + offset + chunk_size),
                              (ghost_array.data() + my_ghosts->second -
                               std::max(ghost_array.data() + my_ghosts->first,
                                        ghost_array_ptr + offset +
                                          chunk_size))),
                            0);
                        }
                    }
                  offset += chunk_size;
                }
              AssertDimension(offset, ghost_targets_data[i].second);
            }

          AssertThrow(
            static_cast<std::size_t>(ghost_targets_data[i].second) *
                sizeof(Number) <
              static_cast<std::size_t>(std::numeric_limits<int>::max()),
            ExcMessage("Index overflow: Maximum message size in MPI is 2GB. "
                       "The number of ghost entries times the size of 'Number' "
                       "exceeds this value. This is not supported."));
          if (std::is_same_v<MemorySpaceType, MemorySpace::Default>)
            Kokkos::fence();
          const int ierr =
            MPI_Isend(ghost_array_ptr,
                      ghost_targets_data[i].second * sizeof(Number),
                      MPI_BYTE,
                      ghost_targets_data[i].first,
                      mpi_tag,
                      communicator,
                      &requests[n_import_targets + i]);
          AssertThrowMPI(ierr);

          ghost_array_ptr += ghost_targets_data[i].second;
        }
    }



    namespace internal
    {
      // In the import_from_ghosted_array_finish we need to invoke abs() also
      // on unsigned data types, which is ill-formed on newer C++
      // standards. To avoid this, we use std::abs on default types but
      // simply return the number on unsigned types
      template <typename Number>
      std::enable_if_t<!std::is_unsigned_v<Number>,
                       typename numbers::NumberTraits<Number>::real_type>
      get_abs(const Number a)
      {
        return std::abs(a);
      }

      template <typename Number>
      std::enable_if_t<std::is_unsigned_v<Number>, Number>
      get_abs(const Number a)
      {
        return a;
      }

      // In the import_from_ghosted_array_finish we might need to calculate the
      // maximal and minimal value for the given number type, which is not
      // straight forward for complex numbers. Therefore, comparison of complex
      // numbers is prohibited and throws an exception.
      template <typename Number>
      DEAL_II_HOST_DEVICE Number
      get_min(const Number a, const Number b)
      {
        return std::min(a, b);
      }

      template <typename Number>
      std::complex<Number>
      get_min(const std::complex<Number> a, const std::complex<Number>)
      {
        AssertThrow(false,
                    ExcMessage("VectorOperation::min not "
                               "implemented for complex numbers"));
        return a;
      }

      template <typename Number>
      DEAL_II_HOST_DEVICE Number
      get_max(const Number a, const Number b)
      {
        return std::max(a, b);
      }

      template <typename Number>
      std::complex<Number>
      get_max(const std::complex<Number> a, const std::complex<Number>)
      {
        AssertThrow(false,
                    ExcMessage("VectorOperation::max not "
                               "implemented for complex numbers"));
        return a;
      }
    } // namespace internal



    template <typename Number, typename MemorySpaceType>
    void
    Partitioner::import_from_ghosted_array_finish(
      const VectorOperation::values                   vector_operation,
      const ArrayView<const Number, MemorySpaceType> &temporary_storage,
      const ArrayView<Number, MemorySpaceType>       &locally_owned_array,
      const ArrayView<Number, MemorySpaceType>       &ghost_array,
      std::vector<MPI_Request>                       &requests) const
    {
      AssertDimension(temporary_storage.size(), n_import_indices());
      Assert(ghost_array.size() == n_ghost_indices() ||
               ghost_array.size() == n_ghost_indices_in_larger_set,
             ExcGhostIndexArrayHasWrongSize(ghost_array.size(),
                                            n_ghost_indices(),
                                            n_ghost_indices_in_larger_set));

      // in optimized mode, no communication was started, so leave the
      // function directly (and only clear ghosts)
      if constexpr (running_in_debug_mode() == false)
        if (vector_operation == VectorOperation::insert)
          {
            Assert(requests.empty(),
                   ExcInternalError(
                     "Did not expect a non-empty communication "
                     "request when inserting. Check that the same "
                     "vector_operation argument was passed to "
                     "import_from_ghosted_array_start as is passed "
                     "to import_from_ghosted_array_finish."));

            Kokkos::deep_copy(
              Kokkos::View<Number *, typename MemorySpaceType::kokkos_space>(
                ghost_array.data(), ghost_array.size()),
              0);
            return;
          }

      // nothing to do when we neither have import nor ghost indices.
      if (n_ghost_indices() == 0 && n_import_indices() == 0)
        return;

      const unsigned int n_import_targets = import_targets_data.size();
      const unsigned int n_ghost_targets  = ghost_targets_data.size();

#    if defined(DEAL_II_MPI_WITH_DEVICE_SUPPORT)
      // When using device-aware MPI, the set of local indices that are ghosts
      // indices on other processors is expanded in arrays. This is for
      // performance reasons as this can significantly decrease the number of
      // kernel launched. The indices are expanded the first time the function
      // is called.
      if ((std::is_same_v<MemorySpaceType, MemorySpace::Default>)&&(
            import_indices_plain_dev.empty()))
        initialize_import_indices_plain_dev();
#    endif

      if (vector_operation != VectorOperation::insert)
        AssertDimension(n_ghost_targets + n_import_targets, requests.size());
      // first wait for the receive to complete
      if (requests.size() > 0 && n_import_targets > 0)
        {
          AssertDimension(locally_owned_array.size(), locally_owned_size());
          const int ierr =
            MPI_Waitall(n_import_targets, requests.data(), MPI_STATUSES_IGNORE);
          AssertThrowMPI(ierr);

          const Number *read_position = temporary_storage.data();
#    if defined(DEAL_II_MPI_WITH_DEVICE_SUPPORT)
          if constexpr (std::is_same_v<MemorySpaceType, MemorySpace::Default>)
            {
              if (vector_operation == VectorOperation::add)
                {
                  for (const auto &import_indices_plain :
                       import_indices_plain_dev)
                    {
                      const auto chunk_size = import_indices_plain.size();

                      using IndexType = decltype(chunk_size);
                      MemorySpace::Default::kokkos_space::execution_space exec;
                      Kokkos::parallel_for(
                        "dealii::fill locally_owned_array, add",
                        Kokkos::RangePolicy<
                          MemorySpace::Default::kokkos_space::execution_space>(
                          exec, 0, chunk_size),
                        KOKKOS_LAMBDA(IndexType idx) {
                          locally_owned_array
                            .data()[import_indices_plain(idx)] +=
                            read_position[idx];
                        });
                      exec.fence();

                      read_position += chunk_size;
                    }
                }
              else if (vector_operation == VectorOperation::min)
                {
                  for (const auto &import_indices_plain :
                       import_indices_plain_dev)
                    {
                      const auto chunk_size = import_indices_plain.size();

                      using IndexType = decltype(chunk_size);
                      MemorySpace::Default::kokkos_space::execution_space exec;
                      Kokkos::parallel_for(
                        "dealii::fill locally_owned_array, min",
                        Kokkos::RangePolicy<
                          MemorySpace::Default::kokkos_space::execution_space>(
                          exec, 0, chunk_size),
                        KOKKOS_LAMBDA(IndexType idx) {
                          locally_owned_array
                            .data()[import_indices_plain(idx)] =
                            internal::get_min(
                              locally_owned_array
                                .data()[import_indices_plain(idx)],
                              read_position[idx]);
                        });
                      exec.fence();

                      read_position += chunk_size;
                    }
                }
              else if (vector_operation == VectorOperation::max)
                {
                  for (const auto &import_indices_plain :
                       import_indices_plain_dev)
                    {
                      const auto chunk_size = import_indices_plain.size();

                      using IndexType = decltype(chunk_size);
                      MemorySpace::Default::kokkos_space::execution_space exec;
                      Kokkos::parallel_for(
                        "dealii::fill locally_owned_array, max",
                        Kokkos::RangePolicy<
                          MemorySpace::Default::kokkos_space::execution_space>(
                          exec, 0, chunk_size),
                        KOKKOS_LAMBDA(IndexType idx) {
                          locally_owned_array
                            .data()[import_indices_plain(idx)] =
                            internal::get_max(
                              locally_owned_array
                                .data()[import_indices_plain(idx)],
                              read_position[idx]);
                        });
                      exec.fence();

                      read_position += chunk_size;
                    }
                }
              else
                {
                  for (const auto &import_indices_plain :
                       import_indices_plain_dev)
                    {
                      // We can't easily assert here, so we just move the
                      // pointer matching the host code.
                      const auto chunk_size = import_indices_plain.size();
                      read_position += chunk_size;
                    }
                }
            }
          else
#    endif
            {
              // If the operation is no insertion, add the imported data to the
              // local values. For insert, nothing is done here (but in debug
              // mode we assert that the specified value is either zero or
              // matches with the ones already present
              if (vector_operation == VectorOperation::add)
                for (const auto &import_range : import_indices_data)
                  for (unsigned int j = import_range.first;
                       j < import_range.second;
                       j++)
                    locally_owned_array[j] += *read_position++;
              else if (vector_operation == VectorOperation::min)
                for (const auto &import_range : import_indices_data)
                  for (unsigned int j = import_range.first;
                       j < import_range.second;
                       j++)
                    {
                      locally_owned_array[j] =
                        internal::get_min(*read_position,
                                          locally_owned_array[j]);
                      ++read_position;
                    }
              else if (vector_operation == VectorOperation::max)
                for (const auto &import_range : import_indices_data)
                  for (unsigned int j = import_range.first;
                       j < import_range.second;
                       j++)
                    {
                      locally_owned_array[j] =
                        internal::get_max(*read_position,
                                          locally_owned_array[j]);
                      ++read_position;
                    }
              else
                for (const auto &import_range : import_indices_data)
                  for (unsigned int j = import_range.first;
                       j < import_range.second;
                       j++, read_position++)
                    // Below we use relatively large precision in units in the
                    // last place (ULP) as this Assert can be easily triggered
                    // in p::d::SolutionTransfer. The rationale is that during
                    // interpolation on two elements sharing the face, values on
                    // this face obtained from each side might be different due
                    // to additions being done in different order. If the local
                    // value is zero, it indicates that the local process has
                    // not set the value during the cell loop and its value can
                    // be safely overridden.
                    Assert(
                      *read_position == Number() ||
                        internal::get_abs(locally_owned_array[j] -
                                          *read_position) <=
                          internal::get_abs(locally_owned_array[j] +
                                            *read_position) *
                            100000. *
                            std::numeric_limits<typename numbers::NumberTraits<
                              Number>::real_type>::epsilon(),
                      typename dealii::LinearAlgebra::distributed::Vector<
                        Number>::ExcNonMatchingElements(*read_position,
                                                        locally_owned_array[j],
                                                        my_pid));
            }

          AssertDimension(read_position - temporary_storage.data(),
                          n_import_indices());
        }

      // wait for the send operations to complete
      if (requests.size() > 0 && n_ghost_targets > 0)
        {
          const int ierr = MPI_Waitall(n_ghost_targets,
                                       &requests[n_import_targets],
                                       MPI_STATUSES_IGNORE);
          AssertThrowMPI(ierr);
        }
      else
        AssertDimension(n_ghost_indices(), 0);

      // clear the ghost array in case we did not yet do that in the _start
      // function
      if (ghost_array.size() > 0)
        {
          Assert(ghost_array.begin() != nullptr, ExcInternalError());

#    if defined(DEAL_II_MPI_WITH_DEVICE_SUPPORT)
          if constexpr (std::is_same_v<MemorySpaceType, MemorySpace::Default>)
            {
              Kokkos::deep_copy(
                Kokkos::View<Number *, MemorySpace::Default::kokkos_space>(
                  ghost_array.data(), n_ghost_indices()),
                Number(0));
            }
          else
#    endif
            {
              std::fill(ghost_array.data(),
                        ghost_array.data() + n_ghost_indices(),
                        Number(0));
            }
        }

      // clear the compress requests
      requests.resize(0);
    }


#  endif // ifdef DEAL_II_WITH_MPI
#endif   // ifndef DOXYGEN

  } // end of namespace MPI

} // namespace Utilities


DEAL_II_NAMESPACE_CLOSE

#endif
