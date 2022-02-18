// ---------------------------------------------------------------------
//
// Copyright (C) 2019 - 2021 by the deal.II authors
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

#include <deal.II/base/config.h>

#include <deal.II/base/mpi.h>
#include <deal.II/base/mpi_compute_index_owner_internal.h>

DEAL_II_NAMESPACE_OPEN

namespace Utilities
{
  namespace MPI
  {
    namespace internal
    {
      /**
       * An internal namespace used for Utilities::MPI::compute_index_owner()
       * and for Utilities::MPI::Partitioner::set_ghost_indices().
       */
      namespace ComputeIndexOwner
      {
        DictionaryPayLoad::DictionaryPayLoad(
          const std::map<unsigned int,
                         std::vector<std::pair<types::global_dof_index,
                                               types::global_dof_index>>>
            &                        buffers,
          std::vector<unsigned int> &actually_owning_ranks,
          const std::pair<types::global_dof_index, types::global_dof_index>
            &                        local_range,
          std::vector<unsigned int> &actually_owning_rank_list)
          : buffers(buffers)
          , actually_owning_ranks(actually_owning_ranks)
          , local_range(local_range)
          , actually_owning_rank_list(actually_owning_rank_list)
        {
          Assert(local_range.first <= local_range.second, ExcInternalError());
        }



        std::vector<unsigned int>
        DictionaryPayLoad::compute_targets()
        {
          std::vector<unsigned int> targets;
          for (const auto &rank_pair : buffers)
            targets.push_back(rank_pair.first);

          return targets;
        }



        void
        DictionaryPayLoad::create_request(
          const unsigned int                               other_rank,
          std::vector<std::pair<types::global_dof_index,
                                types::global_dof_index>> &send_buffer)
        {
          send_buffer = this->buffers.at(other_rank);
        }



        void
        DictionaryPayLoad::answer_request(
          const unsigned int                                     other_rank,
          const std::vector<std::pair<types::global_dof_index,
                                      types::global_dof_index>> &buffer_recv,
          std::vector<unsigned int> &                            request_buffer)
        {
          (void)request_buffer; // not needed


          // process message: loop over all intervals
          for (auto interval : buffer_recv)
            {
#ifdef DEBUG
              for (types::global_dof_index i = interval.first;
                   i < interval.second;
                   i++)
                Assert(
                  actually_owning_ranks[i - local_range.first] ==
                    numbers::invalid_unsigned_int,
                  ExcMessage(
                    "Multiple processes seem to own the same global index. "
                    "A possible reason is that the sets of locally owned "
                    "indices are not distinct."));
              Assert(interval.first < interval.second, ExcInternalError());
              Assert(
                local_range.first <= interval.first &&
                  interval.second <= local_range.second,
                ExcMessage(
                  "The specified interval is not handled by the current process."));
#endif
              std::fill(actually_owning_ranks.data() + interval.first -
                          local_range.first,
                        actually_owning_ranks.data() + interval.second -
                          local_range.first,
                        other_rank);
            }
          actually_owning_rank_list.push_back(other_rank);
        }
      } // namespace ComputeIndexOwner
    }   // namespace internal
  }     // namespace MPI
} // namespace Utilities

DEAL_II_NAMESPACE_CLOSE
