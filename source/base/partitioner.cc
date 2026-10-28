// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2012 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/partitioner.h>
#include <deal.II/base/partitioner.templates.h>

#include <boost/serialization/utility.hpp>

#include <Kokkos_Core.hpp>

#include <algorithm>
#include <limits>
#include <string>
#include <utility>
#include <vector>

DEAL_II_NAMESPACE_OPEN

#ifndef DOXYGEN
namespace Utilities
{
  namespace MPI
  {
    Partitioner::Partitioner()
      : global_size(0)
      , local_range_data(
          std::pair<types::global_dof_index, types::global_dof_index>(0, 0))
      , n_ghost_indices_data(0)
      , n_import_indices_data(0)
      , n_ghost_indices_in_larger_set(0)
      , my_pid(0)
      , n_procs(1)
      , communicator(MPI_COMM_SELF)
      , have_ghost_indices(false)
    {}



    Partitioner::Partitioner(const unsigned int size)
      : global_size(size)
      , locally_owned_range_data(size)
      , local_range_data(
          std::pair<types::global_dof_index, types::global_dof_index>(0, size))
      , n_ghost_indices_data(0)
      , n_import_indices_data(0)
      , n_ghost_indices_in_larger_set(0)
      , my_pid(0)
      , n_procs(1)
      , communicator(MPI_COMM_SELF)
      , have_ghost_indices(false)
    {
      locally_owned_range_data.add_range(0, size);
      locally_owned_range_data.compress();
      ghost_indices_data.set_size(size);
    }



    Partitioner::Partitioner(const types::global_dof_index local_size,
                             const types::global_dof_index ghost_size,
                             const MPI_Comm                communicator)
      : global_size(Utilities::MPI::sum<types::global_dof_index>(local_size,
                                                                 communicator))
      , locally_owned_range_data(global_size)
      , local_range_data{0, local_size}
      , n_ghost_indices_data(ghost_size)
      , n_import_indices_data(0)
      , n_ghost_indices_in_larger_set(0)
      , my_pid(Utilities::MPI::this_mpi_process(communicator))
      , n_procs(Utilities::MPI::n_mpi_processes(communicator))
      , communicator(communicator)
      , have_ghost_indices(true)
    {
      types::global_dof_index prefix_sum = 0;

#  ifdef DEAL_II_WITH_MPI
      const int ierr =
        MPI_Exscan(&local_size,
                   &prefix_sum,
                   1,
                   Utilities::MPI::mpi_type_id_for_type<decltype(prefix_sum)>,
                   MPI_SUM,
                   communicator);
      AssertThrowMPI(ierr);
#  endif

      local_range_data = {prefix_sum, prefix_sum + local_size};

      locally_owned_range_data.add_range(prefix_sum, prefix_sum + local_size);
      locally_owned_range_data.compress();
    }



    Partitioner::Partitioner(const IndexSet &locally_owned_indices,
                             const IndexSet &ghost_indices_in,
                             const MPI_Comm  communicator_in)
      : global_size(
          static_cast<types::global_dof_index>(locally_owned_indices.size()))
      , n_ghost_indices_data(0)
      , n_import_indices_data(0)
      , n_ghost_indices_in_larger_set(0)
      , my_pid(0)
      , n_procs(1)
      , communicator(communicator_in)
      , have_ghost_indices(false)
    {
      set_owned_indices(locally_owned_indices);
      set_ghost_indices(ghost_indices_in);
    }



    Partitioner::Partitioner(const IndexSet &locally_owned_indices,
                             const MPI_Comm  communicator_in)
      : global_size(
          static_cast<types::global_dof_index>(locally_owned_indices.size()))
      , n_ghost_indices_data(0)
      , n_import_indices_data(0)
      , n_ghost_indices_in_larger_set(0)
      , my_pid(0)
      , n_procs(1)
      , communicator(communicator_in)
      , have_ghost_indices(false)
    {
      set_owned_indices(locally_owned_indices);
    }



    void
    Partitioner::reinit(const IndexSet &locally_owned_indices,
                        const IndexSet &ghost_indices,
                        const MPI_Comm  communicator_in)
    {
      have_ghost_indices = false;
      communicator       = communicator_in;
      set_owned_indices(locally_owned_indices);
      set_ghost_indices(ghost_indices);
    }



    void
    Partitioner::set_owned_indices(const IndexSet &locally_owned_indices)
    {
      my_pid  = Utilities::MPI::this_mpi_process(communicator);
      n_procs = Utilities::MPI::n_mpi_processes(communicator);

      // set the local range
      Assert(locally_owned_indices.is_contiguous() == true,
             ExcMessage("The index set specified in locally_owned_indices "
                        "is not contiguous."));
      locally_owned_indices.compress();
      if (locally_owned_indices.n_elements() > 0)
        local_range_data =
          std::pair<types::global_dof_index, types::global_dof_index>(
            locally_owned_indices.nth_index_in_set(0),
            locally_owned_indices.nth_index_in_set(0) +
              locally_owned_indices.n_elements());
      AssertThrow(
        local_range_data.second - local_range_data.first <
          static_cast<types::global_dof_index>(
            std::numeric_limits<unsigned int>::max()),
        ExcMessage(
          "Index overflow: This class supports at most 2^32-1 locally owned vector entries"));
      locally_owned_range_data.set_size(locally_owned_indices.size());
      locally_owned_range_data.add_range(local_range_data.first,
                                         local_range_data.second);
      locally_owned_range_data.compress();

      ghost_indices_data.set_size(locally_owned_indices.size());
    }



    void
    Partitioner::set_ghost_indices(const IndexSet &ghost_indices_in,
                                   const IndexSet &larger_ghost_index_set)
    {
      // Set ghost indices from input. To be sure that no entries from the
      // locally owned range are present, subtract the locally owned indices
      // in any case.
      Assert(ghost_indices_in.n_elements() == 0 ||
               ghost_indices_in.size() == locally_owned_range_data.size(),
             ExcDimensionMismatch(ghost_indices_in.size(),
                                  locally_owned_range_data.size()));

      // The expanded import indices for device-aware MPI are set up lazily
      // on first use. Discard them as the communication pattern changes.
      import_indices_plain_dev   = {};
      import_entries_dev         = {};
      import_entry_offsets_dev   = {};
      import_entry_positions_dev = {};

      ghost_indices_data = ghost_indices_in;
      if (ghost_indices_data.size() != locally_owned_range_data.size())
        ghost_indices_data.set_size(locally_owned_range_data.size());
      ghost_indices_data.subtract_set(locally_owned_range_data);
      ghost_indices_data.compress();
      AssertThrow(
        ghost_indices_data.n_elements() <
          static_cast<types::global_dof_index>(
            std::numeric_limits<unsigned int>::max()),
        ExcMessage(
          "Index overflow: This class supports at most 2^32-1 ghost elements"));
      n_ghost_indices_data = ghost_indices_data.n_elements();

      have_ghost_indices =
        Utilities::MPI::logical_or(n_ghost_indices_data > 0, communicator);

      // In the rest of this function, we determine the point-to-point
      // communication pattern of the partitioner. We make up a list with both
      // the processors the ghost indices actually belong to, and the indices
      // that are locally held but ghost indices of other processors. This
      // allows then to import and export data very easily.

      // find out the end index for each processor and communicate it (this
      // implies the start index for the next processor)
#  ifdef DEAL_II_WITH_MPI
      if (n_procs < 2)
        {
          Assert(ghost_indices_data.n_elements() == 0, ExcInternalError());
          Assert(n_import_indices_data == 0, ExcInternalError());
          Assert(n_ghost_indices_data == 0, ExcInternalError());
          return;
        }

      types::global_dof_index my_size = locally_owned_size();

      // Allow non-zero start index for the vector. Part 1:
      // Assume for now that the index set of rank 0 starts with 0
      // and therefore has an increased size.
      if (my_pid == 0)
        my_size += local_range_data.first;

      types::global_dof_index my_shift = 0;
      {
        const int ierr = MPI_Exscan(
          &my_size,
          &my_shift,
          1,
          Utilities::MPI::mpi_type_id_for_type<types::global_dof_index>,
          MPI_SUM,
          communicator);
        AssertThrowMPI(ierr);
      }

      // Allow non-zero start index for the vector. Part 2:
      // We correct the assumption made above and let the
      // index set of rank 0 actually start from the
      // correct value, i.e. we correct the shift to
      // its start.
      if (my_pid == 0)
        my_shift = local_range_data.first;

      // Fix the index start in case the index set could not give us that
      // information.
      if (local_range_data.first == 0 && my_shift != 0)
        {
          const types::global_dof_index old_locally_owned_size =
            locally_owned_size();
          local_range_data.first  = my_shift;
          local_range_data.second = my_shift + old_locally_owned_size;
        }

      const auto [owning_ranks_of_ghosts, import_data] =
        Utilities::MPI::compute_index_owner_and_requesters(
          locally_owned_range_data, ghost_indices_data, communicator);


      {
        ghost_targets_data = {};

        if (owning_ranks_of_ghosts.size() > 0)
          {
            ghost_targets_data.emplace_back(owning_ranks_of_ghosts[0], 0);
            for (auto i : owning_ranks_of_ghosts)
              {
                Assert(i >= ghost_targets_data.back().first,
                       ExcInternalError(
                         "Expect result of ConsensusAlgorithms::Process to be "
                         "sorted."));
                if (i == ghost_targets_data.back().first)
                  ghost_targets_data.back().second++;
                else
                  ghost_targets_data.emplace_back(i, 1);
              }
          }
      }

      // count import requests and set up the compressed indices
      n_import_indices_data = 0;
      import_targets_data   = {};
      import_targets_data.reserve(import_data.size());
      import_indices_chunks_by_rank_data = {};
      import_indices_chunks_by_rank_data.reserve(import_data.size());
      import_indices_chunks_by_rank_data.resize(1);
      for (const auto &i : import_data)
        if (i.second.n_elements() > 0)
          {
            import_targets_data.emplace_back(i.first, i.second.n_elements());
            n_import_indices_data += i.second.n_elements();
            import_indices_chunks_by_rank_data.push_back(
              import_indices_chunks_by_rank_data.back() +
              i.second.n_intervals());
          }

      // transform import indices to local index space
      import_indices_data = {};
      import_indices_data.reserve(import_indices_chunks_by_rank_data.back());
      for (const auto &i : import_data)
        {
          Assert((i.second & locally_owned_range_data) == i.second,
                 ExcInternalError("Requested indices must be in local range"));
          for (auto interval = i.second.begin_intervals();
               interval != i.second.end_intervals();
               ++interval)
            import_indices_data.emplace_back(*interval->begin() -
                                               local_range_data.first,
                                             interval->last() + 1 -
                                               local_range_data.first);
        }

      if constexpr (running_in_debug_mode())
        {
          // simple check: the number of processors to which we want to send
          // ghosts and the processors to which ghosts reference should be the
          // same
          AssertDimension(
            Utilities::MPI::sum(import_targets_data.size(), communicator),
            Utilities::MPI::sum(ghost_targets_data.size(), communicator));

          // simple check: the number of indices to exchange should match from
          // the ghost indices side and the import indices side
          AssertDimension(
            Utilities::MPI::sum(n_import_indices_data, communicator),
            Utilities::MPI::sum(n_ghost_indices_data, communicator));

          // expensive check that the communication channel is sane -> do a
          // ghost exchange step and see whether the ghost indices sent to us by
          // other processes (ghost_indices) are the same as we hold locally
          // (ghost_indices_ref).
          const std::vector<types::global_dof_index> ghost_indices_ref =
            ghost_indices_data.get_index_vector();
          AssertDimension(ghost_indices_ref.size(), n_ghost_indices());
          std::vector<types::global_dof_index> indices_to_send(
            n_import_indices());
          std::vector<types::global_dof_index> ghost_indices(n_ghost_indices());

          const std::vector<types::global_dof_index> my_indices =
            locally_owned_range_data.get_index_vector();
          std::vector<MPI_Request> requests;
          n_ghost_indices_in_larger_set = n_ghost_indices_data;
          export_to_ghosted_array_start(
            127,
            ArrayView<const types::global_dof_index>(my_indices.data(),
                                                     my_indices.size()),
            make_array_view(indices_to_send),
            make_array_view(ghost_indices),
            requests);
          export_to_ghosted_array_finish(make_array_view(ghost_indices),
                                         requests);
          int       flag = 0;
          const int ierr = MPI_Testall(requests.size(),
                                       requests.data(),
                                       &flag,
                                       MPI_STATUSES_IGNORE);
          AssertThrowMPI(ierr);
          Assert(flag == 1,
                 ExcMessage(
                   "MPI found unfinished requests. Check communication setup"));

          for (unsigned int i = 0; i < ghost_indices.size(); ++i)
            AssertDimension(ghost_indices[i], ghost_indices_ref[i]);
        }

#  endif // #ifdef DEAL_II_WITH_MPI

      if (larger_ghost_index_set.size() == 0)
        {
          ghost_indices_subset_chunks_by_rank_data.clear();
          ghost_indices_subset_data.emplace_back(0, n_ghost_indices());
          n_ghost_indices_in_larger_set = n_ghost_indices_data;
        }
      else
        {
          AssertDimension(larger_ghost_index_set.size(),
                          ghost_indices_data.size());
          Assert(
            (larger_ghost_index_set & locally_owned_range_data).n_elements() ==
              0,
            ExcMessage("Ghost index set should not overlap with owned set."));
          Assert((larger_ghost_index_set & ghost_indices_data) ==
                   ghost_indices_data,
                 ExcMessage("Larger ghost index set must contain the tight "
                            "ghost index set."));

          n_ghost_indices_in_larger_set = larger_ghost_index_set.n_elements();

          // first translate tight ghost indices into indices within the large
          // set:
          std::vector<unsigned int> expanded_numbering;
          for (const dealii::IndexSet::size_type index : ghost_indices_data)
            {
              Assert(larger_ghost_index_set.is_element(index),
                     ExcMessage("The given larger ghost index set must contain "
                                "all indices in the actual index set."));
              Assert(
                larger_ghost_index_set.index_within_set(index) <
                  static_cast<types::global_dof_index>(
                    std::numeric_limits<unsigned int>::max()),
                ExcMessage(
                  "Index overflow: This class supports at most 2^32-1 ghost elements"));
              expanded_numbering.push_back(
                larger_ghost_index_set.index_within_set(index));
            }

          // now rework expanded_numbering into ranges and store in:
          std::vector<std::pair<unsigned int, unsigned int>>
            ghost_indices_subset;
          ghost_indices_subset_chunks_by_rank_data.resize(
            ghost_targets_data.size() + 1);
          // also populate ghost_indices_subset_chunks_by_rank_data
          ghost_indices_subset_chunks_by_rank_data[0] = 0;
          unsigned int shift                          = 0;
          for (unsigned int p = 0; p < ghost_targets_data.size(); ++p)
            {
              unsigned int last_index = numbers::invalid_unsigned_int - 1;
              for (unsigned int ii = 0; ii < ghost_targets_data[p].second; ++ii)
                {
                  const unsigned int i = shift + ii;
                  if (expanded_numbering[i] == last_index + 1)
                    // if contiguous, increment the end of last range:
                    ghost_indices_subset.back().second++;
                  else
                    // otherwise start a new range
                    ghost_indices_subset.emplace_back(expanded_numbering[i],
                                                      expanded_numbering[i] +
                                                        1);
                  last_index = expanded_numbering[i];
                }
              shift += ghost_targets_data[p].second;
              ghost_indices_subset_chunks_by_rank_data[p + 1] =
                ghost_indices_subset.size();
            }
          ghost_indices_subset_data = std::move(ghost_indices_subset);
        }
    }



    bool
    Partitioner::is_compatible(const Partitioner &part) const
    {
      // if the partitioner points to the same memory location as the calling
      // processor
      if (&part == this)
        return true;
#  ifdef DEAL_II_WITH_MPI
      if (Utilities::MPI::job_supports_mpi())
        {
          int       communicators_same = 0;
          const int ierr               = MPI_Comm_compare(part.communicator,
                                            communicator,
                                            &communicators_same);
          AssertThrowMPI(ierr);
          if (!(communicators_same == MPI_IDENT ||
                communicators_same == MPI_CONGRUENT))
            return false;
        }
#  endif
      return (global_size == part.global_size &&
              local_range_data == part.local_range_data &&
              ghost_indices_data == part.ghost_indices_data);
    }



    bool
    Partitioner::is_globally_compatible(const Partitioner &part) const
    {
      return Utilities::MPI::logical_and(is_compatible(part), communicator);
    }



    std::size_t
    Partitioner::memory_consumption() const
    {
      std::size_t memory = MemoryConsumption::memory_consumption(global_size);
      memory += locally_owned_range_data.memory_consumption();
      memory += MemoryConsumption::memory_consumption(local_range_data);
      memory += ghost_indices_data.memory_consumption();
      memory += sizeof(n_ghost_indices_data);
      memory += MemoryConsumption::memory_consumption(ghost_targets_data);
      memory += MemoryConsumption::memory_consumption(import_indices_data);
      memory += sizeof(import_indices_plain_dev) +
                sizeof(unsigned int) * import_indices_plain_dev.size();
      memory += sizeof(import_entries_dev) +
                sizeof(unsigned int) * import_entries_dev.size();
      memory += sizeof(import_entry_offsets_dev) +
                sizeof(unsigned int) * import_entry_offsets_dev.size();
      memory += sizeof(import_entry_positions_dev) +
                sizeof(unsigned int) * import_entry_positions_dev.size();
      memory += MemoryConsumption::memory_consumption(n_import_indices_data);
      memory += MemoryConsumption::memory_consumption(import_targets_data);
      memory += MemoryConsumption::memory_consumption(
        import_indices_chunks_by_rank_data);
      memory +=
        MemoryConsumption::memory_consumption(n_ghost_indices_in_larger_set);
      memory += MemoryConsumption::memory_consumption(
        ghost_indices_subset_chunks_by_rank_data);
      memory +=
        MemoryConsumption::memory_consumption(ghost_indices_subset_data);
      memory += MemoryConsumption::memory_consumption(my_pid);
      memory += MemoryConsumption::memory_consumption(n_procs);
      memory += MemoryConsumption::memory_consumption(communicator);
      memory += MemoryConsumption::memory_consumption(have_ghost_indices);
      return memory;
    }



    void
    Partitioner::initialize_import_indices_plain_dev() const
    {
      // Expand the indices of all import targets on the host. The ranges in
      // import_indices_data are sorted by import target, so the expanded
      // indices are in the same order as the data in the temporary storage
      // used for communication.
      std::vector<unsigned int> import_indices_plain_host;
      import_indices_plain_host.reserve(n_import_indices_data);
      for (const auto &import_range : import_indices_data)
        for (unsigned int j = import_range.first; j < import_range.second; ++j)
          import_indices_plain_host.push_back(j);
      AssertDimension(import_indices_plain_host.size(), n_import_indices_data);

      // An index can be a ghost on several processors. Sort the positions
      // by index, keeping the positions of the same index in increasing
      // order, i.e., in the order of the import targets.
      std::vector<std::pair<unsigned int, unsigned int>> index_and_position;
      index_and_position.reserve(import_indices_plain_host.size());
      for (unsigned int p = 0; p < import_indices_plain_host.size(); ++p)
        index_and_position.emplace_back(import_indices_plain_host[p], p);
      std::sort(index_and_position.begin(), index_and_position.end());

      std::vector<unsigned int> import_entries_host;
      std::vector<unsigned int> import_entry_offsets_host;
      std::vector<unsigned int> import_entry_positions_host;
      import_entry_positions_host.reserve(index_and_position.size());
      for (const auto &[index, position] : index_and_position)
        {
          if (import_entries_host.empty() ||
              import_entries_host.back() != index)
            {
              import_entries_host.push_back(index);
              import_entry_offsets_host.push_back(
                import_entry_positions_host.size());
            }
          import_entry_positions_host.push_back(position);
        }
      import_entry_offsets_host.push_back(import_entry_positions_host.size());

      // Move the indices to the device
      const auto copy_to_device = [](const std::vector<unsigned int> &host,
                                     const std::string               &name) {
        Kokkos::View<unsigned int *, MemorySpace::Default::kokkos_space> device(
          Kokkos::view_alloc(name, Kokkos::WithoutInitializing), host.size());
        Kokkos::deep_copy(device,
                          Kokkos::View<const unsigned int *, Kokkos::HostSpace>(
                            host.data(), host.size()));
        return device;
      };
      import_indices_plain_dev =
        copy_to_device(import_indices_plain_host, "import_indices_plain_dev");
      import_entries_dev =
        copy_to_device(import_entries_host, "import_entries_dev");
      import_entry_offsets_dev =
        copy_to_device(import_entry_offsets_host, "import_entry_offsets_dev");
      import_entry_positions_dev = copy_to_device(import_entry_positions_host,
                                                  "import_entry_positions_dev");
    }

  } // namespace MPI

} // end of namespace Utilities

#endif

// explicit instantiations from .templates.h file
#include "base/partitioner.inst"

DEAL_II_NAMESPACE_CLOSE
