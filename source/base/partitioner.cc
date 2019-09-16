// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2019 by the deal.II authors
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
#include <deal.II/base/partitioner.h>
#include <deal.II/base/partitioner.templates.h>

DEAL_II_NAMESPACE_OPEN

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
    Partitioner::reinit(const IndexSet &vector_space_vector_index_set,
                        const IndexSet &read_write_vector_index_set,
                        const MPI_Comm &communicator_in)
    {
      have_ghost_indices = false;
      communicator       = communicator_in;
      set_owned_indices(vector_space_vector_index_set);
      set_ghost_indices(read_write_vector_index_set);
    }



    void
    Partitioner::set_owned_indices(const IndexSet &locally_owned_indices)
    {
      if (Utilities::MPI::job_supports_mpi() == true)
        {
          my_pid  = Utilities::MPI::this_mpi_process(communicator);
          n_procs = Utilities::MPI::n_mpi_processes(communicator);
        }
      else
        {
          my_pid  = 0;
          n_procs = 1;
        }

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
        Utilities::MPI::sum(n_ghost_indices_data, communicator) > 0;

      // In the rest of this function, we determine the point-to-point
      // communication pattern of the partitioner. We make up a list with both
      // the processors the ghost indices actually belong to, and the indices
      // that are locally held but ghost indices of other processors. This
      // allows then to import and export data very easily.

      // find out the end index for each processor and communicate it (this
      // implies the start index for the next processor)
#ifdef DEAL_II_WITH_MPI
      if (n_procs < 2)
        {
          Assert(ghost_indices_data.n_elements() == 0, ExcInternalError());
          Assert(n_import_indices_data == 0, ExcInternalError());
          Assert(n_ghost_indices_data == 0, ExcInternalError());
          return;
        }

      types::global_dof_index my_size = local_size();
      // Allow non-zero start index for the vector. send this data to all
      // processors
      if (my_pid == 0)
        my_size += local_range_data.first;

      types::global_dof_index my_shift = 0;
      {
        const int ierr = MPI_Exscan(&my_size,
                                    &my_shift,
                                    1,
                                    DEAL_II_DOF_INDEX_MPI_TYPE,
                                    MPI_SUM,
                                    communicator);
        AssertThrowMPI(ierr);
      }
      if (my_shift != local_range_data.first)
        {
          const types::global_dof_index old_local_size = local_size();
          local_range_data.first                       = my_shift;
          local_range_data.second = my_shift + old_local_size;
        }

      std::vector<unsigned int> owning_ranks_of_ghosts(
        ghost_indices_data.n_elements());

      // set up dictionary
      internal::ComputeIndexOwner::ConsensusAlgorithmPayload process(
        locally_owned_range_data,
        ghost_indices_data,
        communicator,
        owning_ranks_of_ghosts,
        /* track origins of ghosts*/ true);

      // read dictionary by communicating with the process who owns the index
      // in the static partition (i.e. in the dictionary). This process
      // returns the actual owner of the index.
      ConsensusAlgorithmSelector<
        std::pair<types::global_dof_index, types::global_dof_index>,
        unsigned int>
        consensus_algorithm(process, communicator);
      consensus_algorithm.run();

      {
        ghost_targets_data = {};

        if (owning_ranks_of_ghosts.size() > 0)
          {
            ghost_targets_data.emplace_back(owning_ranks_of_ghosts[0], 0);
            for (auto i : owning_ranks_of_ghosts)
              {
                Assert(i >= ghost_targets_data.back().first,
                       ExcInternalError(
                         "Expect result of ConsensusAlgorithmProcess to be "
                         "sorted"));
                if (i == ghost_targets_data.back().first)
                  ghost_targets_data.back().second++;
                else
                  ghost_targets_data.emplace_back(i, 1);
              }
          }
      }

      // find how much the individual processes that want import from me
      std::map<unsigned int, IndexSet> import_data = process.get_requesters();

      // count import requests and setup the compressed indices
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

#  ifdef DEBUG

      // simple check: the number of processors to which we want to send
      // ghosts and the processors to which ghosts reference should be the
      // same
      AssertDimension(
        Utilities::MPI::sum(import_targets_data.size(), communicator),
        Utilities::MPI::sum(ghost_targets_data.size(), communicator));

      // simple check: the number of indices to exchange should match from the
      // ghost indices side and the import indices side
      AssertDimension(Utilities::MPI::sum(n_import_indices_data, communicator),
                      Utilities::MPI::sum(n_ghost_indices_data, communicator));

      // expensive check that the communication channel is sane -> do a ghost
      // exchange step and see whether the ghost indices sent to us by other
      // processes (ghost_indices) are the same as we hold locally
      // (ghost_indices_ref).
      std::vector<types::global_dof_index> ghost_indices_ref;
      ghost_indices_data.fill_index_vector(ghost_indices_ref);
      AssertDimension(ghost_indices_ref.size(), n_ghost_indices());
      std::vector<types::global_dof_index> indices_to_send(n_import_indices());
      std::vector<types::global_dof_index> ghost_indices(n_ghost_indices());
      std::vector<types::global_dof_index> my_indices;
      locally_owned_range_data.fill_index_vector(my_indices);
      std::vector<MPI_Request> requests;
      n_ghost_indices_in_larger_set = n_ghost_indices_data;
      export_to_ghosted_array_start(127,
                                    ArrayView<const types::global_dof_index>(
                                      my_indices.data(), my_indices.size()),
                                    make_array_view(indices_to_send),
                                    make_array_view(ghost_indices),
                                    requests);
      export_to_ghosted_array_finish(make_array_view(ghost_indices), requests);
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

#  endif

#endif // #ifdef DEAL_II_WITH_MPI

      if (larger_ghost_index_set.size() == 0)
        {
          ghost_indices_subset_chunks_by_rank_data.clear();
          ghost_indices_subset_data.emplace_back(local_size(),
                                                 local_size() +
                                                   n_ghost_indices());
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
          for (dealii::IndexSet::size_type index : ghost_indices_data)
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
              for (unsigned int ii = 0; ii < ghost_targets_data[p].second; ii++)
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
          ghost_indices_subset_data = ghost_indices_subset;
        }
    }



    bool
    Partitioner::is_compatible(const Partitioner &part) const
    {
      // if the partitioner points to the same memory location as the calling
      // processor
      if (&part == this)
        return true;
#ifdef DEAL_II_WITH_MPI
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
#endif
      return (global_size == part.global_size &&
              local_range_data == part.local_range_data &&
              ghost_indices_data == part.ghost_indices_data);
    }



    bool
    Partitioner::is_globally_compatible(const Partitioner &part) const
    {
      return Utilities::MPI::min(static_cast<int>(is_compatible(part)),
                                 communicator) == 1;
    }



    std::size_t
    Partitioner::memory_consumption() const
    {
      std::size_t memory = (3 * sizeof(types::global_dof_index) +
                            4 * sizeof(unsigned int) + sizeof(MPI_Comm));
      memory += MemoryConsumption::memory_consumption(locally_owned_range_data);
      memory += MemoryConsumption::memory_consumption(ghost_targets_data);
      memory += MemoryConsumption::memory_consumption(import_targets_data);
      memory += MemoryConsumption::memory_consumption(import_indices_data);
      memory += MemoryConsumption::memory_consumption(
        import_indices_chunks_by_rank_data);
      memory += MemoryConsumption::memory_consumption(
        ghost_indices_subset_chunks_by_rank_data);
      memory +=
        MemoryConsumption::memory_consumption(ghost_indices_subset_data);
      memory += MemoryConsumption::memory_consumption(ghost_indices_data);
      return memory;
    }

  } // namespace MPI

} // end of namespace Utilities



// explicit instantiations from .templates.h file
#include "partitioner.inst"

DEAL_II_NAMESPACE_CLOSE
