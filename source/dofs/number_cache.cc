// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#include <deal.II/base/memory_consumption.h>
#include <deal.II/dofs/number_cache.h>

#include <numeric>


DEAL_II_NAMESPACE_OPEN

namespace internal
{
  namespace DoFHandlerImplementation
  {
    NumberCache::NumberCache ()
      :
      n_global_dofs (0),
      n_locally_owned_dofs (0)
    {}



    NumberCache::NumberCache (const types::global_dof_index n_global_dofs)
      :
      n_global_dofs (n_global_dofs),
      n_locally_owned_dofs (n_global_dofs),
      locally_owned_dofs (complete_index_set(n_global_dofs)),
      n_locally_owned_dofs_per_processor (1, n_global_dofs),
      locally_owned_dofs_per_processor (1, complete_index_set(n_global_dofs))
    {}



    NumberCache::NumberCache (const std::vector<IndexSet> &locally_owned_dofs_per_processor,
                              const unsigned int           my_rank)
      :
      locally_owned_dofs_per_processor (locally_owned_dofs_per_processor)
    {
      const unsigned int n_procs = locally_owned_dofs_per_processor.size();

      // compress IndexSet representation before using it for anything else
      for (unsigned int p=0; p<n_procs; ++p)
        locally_owned_dofs_per_processor[p].compress();

      n_locally_owned_dofs_per_processor.resize (n_procs);
      for (unsigned int p=0; p<n_procs; ++p)
        n_locally_owned_dofs_per_processor[p]
          = locally_owned_dofs_per_processor[p].n_elements();
      n_locally_owned_dofs = n_locally_owned_dofs_per_processor[my_rank];
      locally_owned_dofs = locally_owned_dofs_per_processor[my_rank];


      n_global_dofs = std::accumulate (n_locally_owned_dofs_per_processor.begin(),
                                       n_locally_owned_dofs_per_processor.end(),
                                       types::global_dof_index(0));
    }



    void
    NumberCache::clear ()
    {
      n_global_dofs = 0;
      n_locally_owned_dofs = 0;
      locally_owned_dofs.clear();
      n_locally_owned_dofs_per_processor.clear();
      locally_owned_dofs_per_processor.clear();
    }

    std::size_t
    NumberCache::memory_consumption () const
    {
      return
        MemoryConsumption::memory_consumption (n_global_dofs) +
        MemoryConsumption::memory_consumption (n_locally_owned_dofs) +
        MemoryConsumption::memory_consumption (locally_owned_dofs) +
        MemoryConsumption::memory_consumption (n_locally_owned_dofs_per_processor) +
        MemoryConsumption::memory_consumption (locally_owned_dofs_per_processor);
    }
  }
}

DEAL_II_NAMESPACE_CLOSE
