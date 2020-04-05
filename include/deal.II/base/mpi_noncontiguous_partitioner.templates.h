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

#ifndef dealii_mpi_noncontiguous_vector_template_h
#define dealii_mpi_noncontiguous_vector_template_h

#include <deal.II/base/config.h>

#include <deal.II/base/mpi.h>
#include <deal.II/base/mpi.templates.h>
#include <deal.II/base/mpi_compute_index_owner_internal.h>
#include <deal.II/base/mpi_noncontiguous_partitioner.h>

#include <deal.II/lac/communication_pattern_base.h>
#include <deal.II/lac/vector_space_vector.h>


DEAL_II_NAMESPACE_OPEN

namespace Utilities
{
  namespace MPI
  {
    template <typename Number>
    NoncontiguousPartitioner<Number>::NoncontiguousPartitioner(
      const IndexSet &indexset_has,
      const IndexSet &indexset_want,
      const MPI_Comm &communicator)
    {
      this->reinit(indexset_has, indexset_want, communicator);
    }



    template <typename Number>
    NoncontiguousPartitioner<Number>::NoncontiguousPartitioner(
      const std::vector<types::global_dof_index> &indices_has,
      const std::vector<types::global_dof_index> &indices_want,
      const MPI_Comm &                            communicator)
    {
      this->reinit(indices_has, indices_want, communicator);
    }



    template <typename Number>
    std::pair<unsigned int, unsigned int>
    NoncontiguousPartitioner<Number>::n_targets()
    {
      return {send_ranks.size(), recv_ranks.size()};
    }



    template <typename Number>
    types::global_dof_index
    NoncontiguousPartitioner<Number>::memory_consumption()
    {
      return MemoryConsumption::memory_consumption(send_ranks) +
             MemoryConsumption::memory_consumption(send_ptr) +
             MemoryConsumption::memory_consumption(send_indices) +
             MemoryConsumption::memory_consumption(send_buffers) +
             MemoryConsumption::memory_consumption(send_requests) +
             MemoryConsumption::memory_consumption(recv_ranks) +
             MemoryConsumption::memory_consumption(recv_ptr) +
             MemoryConsumption::memory_consumption(recv_indices) +
             MemoryConsumption::memory_consumption(recv_buffers) +
             MemoryConsumption::memory_consumption(recv_requests);
    }



    template <typename Number>
    const MPI_Comm &
    NoncontiguousPartitioner<Number>::get_mpi_communicator() const
    {
      return communicator;
    }



    template <typename Number>
    void
    NoncontiguousPartitioner<Number>::reinit(const IndexSet &indexset_has,
                                             const IndexSet &indexset_want,
                                             const MPI_Comm &communicator)
    {
      this->communicator = communicator;

      // clean up
      send_ranks.clear();
      send_ptr.clear();
      send_indices.clear();
      send_buffers.clear();
      send_requests.clear();
      recv_ranks.clear();
      recv_ptr.clear();
      recv_indices.clear();
      recv_buffers.clear();
      recv_requests.clear();

      // setup communication pattern
      std::vector<unsigned int> owning_ranks_of_ghosts(
        indexset_want.n_elements());

      // set up dictionary
      Utilities::MPI::internal::ComputeIndexOwner::ConsensusAlgorithmsPayload
        process(indexset_has,
                indexset_want,
                communicator,
                owning_ranks_of_ghosts,
                true);

      Utilities::MPI::ConsensusAlgorithms::Selector<
        std::pair<types::global_dof_index, types::global_dof_index>,
        unsigned int>
        consensus_algorithm(process, communicator);
      consensus_algorithm.run();

      // setup map of processes from where this rank will receive values
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

            for (const auto &cell_index : target_with_indexset.second)
              recv_indices.push_back(cell_index);

            recv_ptr.push_back(recv_indices.size());
          }

        recv_buffers.resize(recv_indices.size());
        recv_requests.resize(recv_map.size());
      }

      {
        const auto targets_with_indexset = process.get_requesters();

        send_ptr.push_back(send_indices.size() /*=0*/);
        for (const auto &target_with_indexset : targets_with_indexset)
          {
            send_ranks.push_back(target_with_indexset.first);

            for (const auto &cell_index : target_with_indexset.second)
              send_indices.push_back(indexset_has.index_within_set(cell_index));

            send_ptr.push_back(send_indices.size());
          }

        send_buffers.resize(send_indices.size());
        send_requests.resize(targets_with_indexset.size());
      }
    }



    template <typename Number>
    void
    NoncontiguousPartitioner<Number>::reinit(
      const std::vector<types::global_dof_index> &indices_has,
      const std::vector<types::global_dof_index> &indices_want,
      const MPI_Comm &                            communicator)
    {
      // step 0) clean vectors from numbers::invalid_dof_index (indicating
      //         padding)
      std::vector<types::global_dof_index> indices_has_clean;
      indices_has_clean.reserve(indices_has.size());

      for (const auto &i : indices_has)
        if (i != numbers::invalid_dof_index)
          indices_has_clean.push_back(i);

      std::vector<types::global_dof_index> indices_want_clean;
      indices_want_clean.reserve(indices_want.size());

      for (const auto &i : indices_want)
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

      // step 2) setup internal data structures with indexset
      this->reinit(index_set_has, index_set_want, communicator);

      // step 3) fix inner data structures so that it is sorted as
      //         in the original vector
      {
        std::vector<types::global_dof_index> temp_map_send(
          index_set_has.n_elements());

        for (types::global_dof_index i = 0; i < indices_has.size(); i++)
          if (indices_has[i] != numbers::invalid_dof_index)
            temp_map_send[index_set_has.index_within_set(indices_has[i])] = i;

        for (auto &i : send_indices)
          i = temp_map_send[i];
      }

      {
        std::vector<types::global_dof_index> temp_map_recv(
          index_set_want.n_elements());

        for (types::global_dof_index i = 0; i < indices_want.size(); i++)
          if (indices_want[i] != numbers::invalid_dof_index)
            temp_map_recv[index_set_want.index_within_set(indices_want[i])] = i;

        for (auto &i : recv_indices)
          i = temp_map_recv[i];
      }
    }


    template <typename Number>
    template <typename VectorType>
    void
    NoncontiguousPartitioner<Number>::update_values(VectorType &      dst,
                                                    const VectorType &src) const
    {
      const auto tag = internal::Tags::noncontiguous_partitioner_update_values;

      this->update_values_start(src, tag);
      this->update_values_finish(dst, tag);
    }



    template <typename Number>
    template <typename VectorType>
    void
    NoncontiguousPartitioner<Number>::update_values_start(
      const VectorType & src,
      const unsigned int tag) const
    {
#ifndef DEAL_II_WITH_MPI
      (void)src;
      (void)tag;
      Assert(false, ExcNeedsMPI());
#else
      // post recv
      for (types::global_dof_index i = 0; i < recv_ranks.size(); i++)
        {
          const auto ierr = MPI_Irecv(recv_buffers.data() + recv_ptr[i],
                                      recv_ptr[i + 1] - recv_ptr[i],
                                      Utilities::MPI::internal::mpi_type_id(
                                        recv_buffers.data()),
                                      recv_ranks[i],
                                      tag,
                                      communicator,
                                      &recv_requests[i]);
          AssertThrowMPI(ierr);
        }

      auto src_iterator = src.begin();

      // post send
      for (types::global_dof_index i = 0; i < send_ranks.size(); i++)
        {
          // collect data to be send
          for (types::global_dof_index j = send_ptr[i], c = 0;
               j < send_ptr[i + 1];
               j++)
            send_buffers[send_ptr[i] + c++] = src_iterator[send_indices[j]];

          // send data
          const auto ierr = MPI_Isend(send_buffers.data() + send_ptr[i],
                                      send_ptr[i + 1] - send_ptr[i],
                                      Utilities::MPI::internal::mpi_type_id(
                                        send_buffers.data()),
                                      send_ranks[i],
                                      tag,
                                      communicator,
                                      &send_requests[i]);
          AssertThrowMPI(ierr);
        }
#endif
    }



    template <typename Number>
    template <typename VectorType>
    void
    NoncontiguousPartitioner<Number>::update_values_finish(
      VectorType &       dst,
      const unsigned int tag) const
    {
      (void)tag;

#ifndef DEAL_II_WITH_MPI
      (void)dst;
      Assert(false, ExcNeedsMPI());
#else
      auto dst_iterator = dst.begin();

      // receive all data packages and copy data from buffers
      for (types::global_dof_index proc = 0; proc < recv_requests.size();
           proc++)
        {
          int        i;
          MPI_Status status;
          const auto ierr = MPI_Waitany(recv_requests.size(),
                                        recv_requests.data(),
                                        &i,
                                        &status);
          AssertThrowMPI(ierr);

          for (types::global_dof_index j = recv_ptr[i], c = 0;
               j < recv_ptr[i + 1];
               j++)
            dst_iterator[recv_indices[j]] = recv_buffers[recv_ptr[i] + c++];
        }

      // wait that all data packages have been sent
      const auto ierr = MPI_Waitall(send_requests.size(),
                                    send_requests.data(),
                                    MPI_STATUSES_IGNORE);
      AssertThrowMPI(ierr);
#endif
    }

  } // namespace MPI
} // namespace Utilities

DEAL_II_NAMESPACE_CLOSE

#endif
