// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2019 by the deal.II authors
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

#include <deal.II/base/memory_consumption.h>

#include <deal.II/grid/tria_levels.h>

DEAL_II_NAMESPACE_OPEN


namespace internal
{
  namespace TriangulationImplementation
  {
    std::size_t
    TriaLevel::memory_consumption() const
    {
      return (MemoryConsumption::memory_consumption(refine_flags) +
              MemoryConsumption::memory_consumption(coarsen_flags) +
              MemoryConsumption::memory_consumption(active_cell_indices) +
              MemoryConsumption::memory_consumption(neighbors) +
              MemoryConsumption::memory_consumption(subdomain_ids) +
              MemoryConsumption::memory_consumption(level_subdomain_ids) +
              MemoryConsumption::memory_consumption(parents) +
              MemoryConsumption::memory_consumption(direction_flags) +
              MemoryConsumption::memory_consumption(cells));
    }
  } // namespace TriangulationImplementation
} // namespace internal

DEAL_II_NAMESPACE_CLOSE
