// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/mpi.h>

#include <deal.II/dofs/number_cache.h>

#include <numeric>


DEAL_II_NAMESPACE_OPEN

namespace internal
{
  namespace DoFHandlerImplementation
  {
    NumberCache::NumberCache()
      : NumberCache(0)
    {}



    NumberCache::NumberCache(const types::global_dof_index n_global_dofs)
      : n_global_dofs(n_global_dofs)
      , n_locally_owned_dofs(n_global_dofs)
      , locally_owned_dofs(complete_index_set(n_global_dofs))
      , n_locally_owned_dofs_per_processor(1, n_global_dofs)
      , locally_owned_dofs_per_processor(1, complete_index_set(n_global_dofs))
    {}



    NumberCache::NumberCache(
      const std::vector<IndexSet> &locally_owned_dofs_per_processor,
      const unsigned int           my_rank)
      : locally_owned_dofs_per_processor(locally_owned_dofs_per_processor)
    {
      const unsigned int n_procs = locally_owned_dofs_per_processor.size();

      // compress IndexSet representation before using it for anything else
      for (unsigned int p = 0; p < n_procs; ++p)
        locally_owned_dofs_per_processor[p].compress();

      n_locally_owned_dofs_per_processor.resize(n_procs);
      for (unsigned int p = 0; p < n_procs; ++p)
        n_locally_owned_dofs_per_processor[p] =
          locally_owned_dofs_per_processor[p].n_elements();
      n_locally_owned_dofs = n_locally_owned_dofs_per_processor[my_rank];
      locally_owned_dofs   = locally_owned_dofs_per_processor[my_rank];


      n_global_dofs =
        std::accumulate(n_locally_owned_dofs_per_processor.begin(),
                        n_locally_owned_dofs_per_processor.end(),
                        types::global_dof_index(0));
    }



    void
    NumberCache::clear()
    {
      n_global_dofs        = 0;
      n_locally_owned_dofs = 0;
      locally_owned_dofs.clear();
      n_locally_owned_dofs_per_processor.clear();
      locally_owned_dofs_per_processor.clear();
    }



    std::vector<types::global_dof_index>
    NumberCache::get_n_locally_owned_dofs_per_processor(
      const MPI_Comm mpi_communicator) const
    {
      if (n_global_dofs == 0)
        return std::vector<types::global_dof_index>();
      else if (n_locally_owned_dofs_per_processor.empty() == false)
        {
          AssertDimension(n_locally_owned_dofs_per_processor.size(),
                          Utilities::MPI::n_mpi_processes(mpi_communicator));
          return n_locally_owned_dofs_per_processor;
        }
      else
        {
          return Utilities::MPI::all_gather(mpi_communicator,
                                            n_locally_owned_dofs);
        }
    }



    std::vector<IndexSet>
    NumberCache::get_locally_owned_dofs_per_processor(
      const MPI_Comm mpi_communicator) const
    {
      AssertDimension(locally_owned_dofs.size(), n_global_dofs);
      if (n_global_dofs == 0)
        return std::vector<IndexSet>();
      else if (locally_owned_dofs_per_processor.empty() == false)
        {
          AssertDimension(locally_owned_dofs_per_processor.size(),
                          Utilities::MPI::n_mpi_processes(mpi_communicator));
          return locally_owned_dofs_per_processor;
        }
      else
        {
          return Utilities::MPI::all_gather(mpi_communicator,
                                            locally_owned_dofs);
        }
    }


    std::size_t
    NumberCache::memory_consumption() const
    {
      return MemoryConsumption::memory_consumption(n_global_dofs) +
             MemoryConsumption::memory_consumption(n_locally_owned_dofs) +
             MemoryConsumption::memory_consumption(locally_owned_dofs) +
             MemoryConsumption::memory_consumption(
               n_locally_owned_dofs_per_processor) +
             MemoryConsumption::memory_consumption(
               locally_owned_dofs_per_processor);
    }
  } // namespace DoFHandlerImplementation
} // namespace internal

DEAL_II_NAMESPACE_CLOSE
