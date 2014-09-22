// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2013 by the deal.II authors
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

DEAL_II_NAMESPACE_OPEN

namespace internal
{
  namespace DoFHandler
  {
    NumberCache::NumberCache ()
      :
      n_global_dofs (0),
      n_locally_owned_dofs (0)
    {}


    void NumberCache::clear ()
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
