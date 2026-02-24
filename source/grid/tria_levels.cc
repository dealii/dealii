// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2006 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

#include <deal.II/base/memory_consumption.h>

#include <deal.II/grid/tria_levels.h>

DEAL_II_NAMESPACE_OPEN


namespace internal
{
  namespace TriangulationImplementation
  {
    template <int dim, int spacedim>
    std::size_t
    TriaLevel<dim, spacedim>::memory_consumption() const
    {
      return (
        MemoryConsumption::memory_consumption(refine_flags) +
        MemoryConsumption::memory_consumption(refine_choice) +
        MemoryConsumption::memory_consumption(coarsen_flags) +
        MemoryConsumption::memory_consumption(active_cell_indices) +
        MemoryConsumption::memory_consumption(global_active_cell_indices) +
        MemoryConsumption::memory_consumption(global_level_cell_indices) +
        MemoryConsumption::memory_consumption(neighbors) +
        MemoryConsumption::memory_consumption(subdomain_ids) +
        MemoryConsumption::memory_consumption(level_subdomain_ids) +
        MemoryConsumption::memory_consumption(parents) +
        MemoryConsumption::memory_consumption(direction_flags) +
        MemoryConsumption::memory_consumption(cells) +
        MemoryConsumption::memory_consumption(face_orientations) +
        MemoryConsumption::memory_consumption(reference_cell) +
        MemoryConsumption::memory_consumption(cell_vertex_indices_cache));
    }

    // explicit instantiations; note: we need them all for all dimensions

    template class TriaLevel<1, 1>;
    template class TriaLevel<1, 2>;
    template class TriaLevel<1, 3>;
    template class TriaLevel<2, 2>;
    template class TriaLevel<2, 3>;
    template class TriaLevel<3, 3>;
  } // namespace TriangulationImplementation
} // namespace internal

DEAL_II_NAMESPACE_CLOSE
