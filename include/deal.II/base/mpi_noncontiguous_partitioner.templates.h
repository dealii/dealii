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

#ifndef dealii_mpi_noncontiguous_partitioner_templates_h
#define dealii_mpi_noncontiguous_partitioner_templates_h

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
    NoncontiguousPartitioner::NoncontiguousPartitioner(
      const IndexSet &indexset_has,
      const IndexSet &indexset_want,
      const MPI_Comm &communicator)
    {
      this->reinit(indexset_has, indexset_want, communicator);
    }



    NoncontiguousPartitioner::NoncontiguousPartitioner(
      const std::vector<types::global_dof_index> &indices_has,
      const std::vector<types::global_dof_index> &indices_want,
      const MPI_Comm &                            communicator)
    {
      this->reinit(indices_has, indices_want, communicator);
    }



    std::pair<unsigned int, unsigned int>
    NoncontiguousPartitioner::n_targets()
    {
      return {send_ranks.size(), recv_ranks.size()};
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



    const MPI_Comm &
    NoncontiguousPartitioner::get_mpi_communicator() const
    {
      return communicator;
    }



    void
    NoncontiguousPartitioner::reinit(const IndexSet &indexset_has,
                                     const IndexSet &indexset_want,
                                     const MPI_Comm &communicator)
    {
      this->communicator = communicator;

      // clean up
      send_ranks.clear();
      send_ptr.clear();
      send_indices.clear();
      recv_ranks.clear();
      recv_ptr.clear();
      recv_indices.clear();
      buffers.clear();
      requests.clear();

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
      }

      {
        const auto targets_with_indexset = process.get_requesters();

        send_ptr.push_back(recv_ptr.back());
        for (const auto &target_with_indexset : targets_with_indexset)
          {
            send_ranks.push_back(target_with_indexset.first);

            for (const auto &cell_index : target_with_indexset.second)
              send_indices.push_back(indexset_has.index_within_set(cell_index));

            send_ptr.push_back(send_indices.size() + recv_ptr.back());
          }
      }
    }



    void
    NoncontiguousPartitioner::reinit(
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
    void
    NoncontiguousPartitioner::export_to_ghosted_array(
      const ArrayView<const Number> &src,
      const ArrayView<Number> &      dst) const
    {
      // allocate internal memory since needed
      if (requests.size() != send_ranks.size() + recv_ranks.size())
        requests.resize(send_ranks.size() + recv_ranks.size());

      if (this->buffers.size() != send_ptr.back() * sizeof(Number))
        this->buffers.resize(send_ptr.back() * sizeof(Number), 0);

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
      AssertIndexRange(communication_channel, 10);

      const auto tag =
        communication_channel +
        internal::Tags::noncontiguous_partitioner_update_ghost_values;

      // post recv
      for (types::global_dof_index i = 0; i < recv_ranks.size(); i++)
        {
          const auto ierr =
            MPI_Irecv(buffers.data() + recv_ptr[i],
                      recv_ptr[i + 1] - recv_ptr[i],
                      Utilities::MPI::internal::mpi_type_id(buffers.data()),
                      recv_ranks[i],
                      tag,
                      communicator,
                      &requests[i + send_ranks.size()]);
          AssertThrowMPI(ierr);
        }

      auto src_iterator = src.begin();

      // post send
      for (types::global_dof_index i = 0, k = 0; i < send_ranks.size(); i++)
        {
          // collect data to be send
          for (types::global_dof_index j = send_ptr[i]; j < send_ptr[i + 1];
               j++)
            buffers[j] = src_iterator[send_indices[k++]];

          // send data
          const auto ierr =
            MPI_Isend(buffers.data() + send_ptr[i],
                      send_ptr[i + 1] - send_ptr[i],
                      Utilities::MPI::internal::mpi_type_id(buffers.data()),
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
      auto dst_iterator = dst.begin();

      // receive all data packages and copy data from buffers
      for (types::global_dof_index proc = 0; proc < recv_ranks.size(); proc++)
        {
          int        i;
          MPI_Status status;
          const auto ierr = MPI_Waitany(recv_ranks.size(),
                                        requests.data() + send_ranks.size(),
                                        &i,
                                        &status);
          AssertThrowMPI(ierr);

          for (types::global_dof_index j = recv_ptr[i], c = 0;
               j < recv_ptr[i + 1];
               j++)
            dst_iterator[recv_indices[j]] = buffers[recv_ptr[i] + c++];
        }

      // wait that all data packages have been sent
      const auto ierr =
        MPI_Waitall(send_ranks.size(), requests.data(), MPI_STATUSES_IGNORE);
      AssertThrowMPI(ierr);
#endif
    }

  } // namespace MPI
} // namespace Utilities

DEAL_II_NAMESPACE_CLOSE

#endif
