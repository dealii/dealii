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

#include <deal.II/base/array_view.h>
#include <deal.II/base/mpi_noncontiguous_partitioner.h>
#include <deal.II/base/mpi_noncontiguous_partitioner.templates.h>

#include <deal.II/lac/la_parallel_vector.h>

#include <boost/serialization/utility.hpp>


DEAL_II_NAMESPACE_OPEN

namespace Utilities
{
  namespace MPI
  {
    NoncontiguousPartitioner::NoncontiguousPartitioner(
      const IndexSet &indexset_has,
      const IndexSet &indexset_want,
      const MPI_Comm  communicator)
    {
      this->reinit(indexset_has, indexset_want, communicator);
    }



    NoncontiguousPartitioner::NoncontiguousPartitioner(
      const std::vector<types::global_dof_index> &indices_has,
      const std::vector<types::global_dof_index> &indices_want,
      const MPI_Comm                              communicator)
    {
      this->reinit(indices_has, indices_want, communicator);
    }



    std::pair<unsigned int, unsigned int>
    NoncontiguousPartitioner::n_targets() const
    {
      return {send_ranks.size(), recv_ranks.size()};
    }



    unsigned int
    NoncontiguousPartitioner::temporary_storage_size() const
    {
      return send_ptr.back();
    }



    types::global_dof_index
    NoncontiguousPartitioner::memory_consumption()
    {
      return MemoryConsumption::memory_consumption(send_ranks) +
             MemoryConsumption::memory_consumption(send_ptr) +
             MemoryConsumption::memory_consumption(send_indices) +
             MemoryConsumption::memory_consumption(recv_ranks) +
             MemoryConsumption::memory_consumption(recv_ptr) +
             MemoryConsumption::memory_consumption(recv_indices) +
             MemoryConsumption::memory_consumption(buffers) +
             MemoryConsumption::memory_consumption(requests);
    }



    MPI_Comm
    NoncontiguousPartitioner::get_mpi_communicator() const
    {
      return communicator;
    }



    void
    NoncontiguousPartitioner::reinit(const IndexSet &indexset_has,
                                     const IndexSet &indexset_want,
                                     const MPI_Comm  communicator)
    {
      this->communicator = communicator;

      // clean up
      send_ranks.clear();
      send_ptr.clear();
      send_indices.clear();
      recv_ranks.clear();
      recv_ptr.clear();
      recv_indices.clear();
      recv_indices_duplicates.clear();
      recv_indices_duplicates_ptr.clear();
      buffers.clear();
      requests.clear();

      // set up communication pattern
      const auto [owning_ranks_of_ghosts, targets_with_indexset] =
        Utilities::MPI::compute_index_owner_and_requesters(indexset_has,
                                                           indexset_want,
                                                           communicator);

      // set up map of processes from where this rank will receive values
      {
        std::map<unsigned int, std::vector<types::global_dof_index>> recv_map;

        for (const auto &owner : owning_ranks_of_ghosts)
          recv_map[owner] = std::vector<types::global_dof_index>();

        for (types::global_dof_index i = 0; i < owning_ranks_of_ghosts.size();
             i++)
          recv_map[owning_ranks_of_ghosts[i]].push_back(i);

        recv_ptr.push_back(recv_indices.size() /*=0*/);
        for (const auto &target_with_indexset : recv_map)
          {
            recv_ranks.push_back(target_with_indexset.first);

            for (const auto cell_index : target_with_indexset.second)
              recv_indices.push_back(cell_index);

            recv_ptr.push_back(recv_indices.size());
          }
      }

      {
        send_ptr.push_back(recv_ptr.back());
        for (const auto &target_with_indexset : targets_with_indexset)
          {
            send_ranks.push_back(target_with_indexset.first);

            for (const auto cell_index : target_with_indexset.second)
              send_indices.push_back(indexset_has.index_within_set(cell_index));

            send_ptr.push_back(send_indices.size() + recv_ptr.back());
          }
      }
    }



    void
    NoncontiguousPartitioner::reinit(
      const std::vector<types::global_dof_index> &indices_has,
      const std::vector<types::global_dof_index> &indices_want,
      const MPI_Comm                              communicator)
    {
      // step 0) clean vectors from numbers::invalid_dof_index (indicating
      //         padding)
      std::vector<types::global_dof_index> indices_has_clean;
      indices_has_clean.reserve(indices_has.size());

      for (const auto i : indices_has)
        if (i != numbers::invalid_dof_index)
          indices_has_clean.push_back(i);

      std::vector<types::global_dof_index> indices_want_clean;
      indices_want_clean.reserve(indices_want.size());

      for (const auto i : indices_want)
        if (i != numbers::invalid_dof_index)
          indices_want_clean.push_back(i);

      // step 0) determine "number of degrees of freedom" needed for IndexSet
      const types::global_dof_index local_n_dofs_has =
        indices_has_clean.empty() ?
          0 :
          (*std::max_element(indices_has_clean.begin(),
                             indices_has_clean.end()) +
           1);

      const types::global_dof_index local_n_dofs_want =
        indices_want_clean.empty() ?
          0 :
          (*std::max_element(indices_want_clean.begin(),
                             indices_want_clean.end()) +
           1);

      const types::global_dof_index n_dofs =
        Utilities::MPI::max(std::max(local_n_dofs_has, local_n_dofs_want),
                            communicator);

      // step 1) convert vectors to indexsets (sorted!)
      IndexSet index_set_has(n_dofs);
      index_set_has.add_indices(indices_has_clean.begin(),
                                indices_has_clean.end());

      IndexSet index_set_want(n_dofs);
      index_set_want.add_indices(indices_want_clean.begin(),
                                 indices_want_clean.end());

      // step 2) set up internal data structures with indexset
      this->reinit(index_set_has, index_set_want, communicator);

      // step 3) fix inner data structures so that it is sorted as
      //         in the original vector
      {
        std::vector<types::global_dof_index> temp_map_send(
          index_set_has.n_elements());

        for (types::global_dof_index i = 0; i < indices_has.size(); ++i)
          if (indices_has[i] != numbers::invalid_dof_index)
            temp_map_send[index_set_has.index_within_set(indices_has[i])] = i;

        for (auto &i : send_indices)
          i = temp_map_send[i];
      }

      {
        std::vector<std::vector<types::global_dof_index>> temp_map_recv(
          index_set_want.n_elements());

        for (types::global_dof_index i = 0; i < indices_want.size(); ++i)
          if (indices_want[i] != numbers::invalid_dof_index)
            temp_map_recv[index_set_want.index_within_set(indices_want[i])]
              .push_back(i);

        const bool use_fast_path =
          std::all_of(temp_map_recv.begin(),
                      temp_map_recv.end(),
                      [&](const auto &x) { return x.size() == 1; });

        if (use_fast_path)
          {
            for (auto &i : recv_indices)
              i = temp_map_recv[i][0];
          }
        else
          {
            recv_indices_duplicates_ptr = {0};

            for (const auto &indices : temp_map_recv)
              {
                for (const auto &index : indices)
                  recv_indices_duplicates.emplace_back(index);

                recv_indices_duplicates_ptr.push_back(
                  recv_indices_duplicates.size());
              }
          }
      }
    }
  } // namespace MPI
} // namespace Utilities

#include "base/mpi_noncontiguous_partitioner.inst"

DEAL_II_NAMESPACE_CLOSE
