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
    TriaLevel<dim, spacedim>::TriaLevel()
      : children_per_object(numbers::invalid_unsigned_int)
      , faces_per_object(numbers::invalid_unsigned_int)
    {}



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



    template <int dim, int spacedim>
    void
    TriaLevel<dim, spacedim>::allocate(const std::size_t n_cells,
                                       const bool        orientation_needed)
    {
      Assert(children_per_object == cells.children_per_object,
             ExcInternalError());
      Assert(faces_per_object == cells.faces_per_object, ExcInternalError());
      active_cell_indices.assign(n_cells, numbers::invalid_unsigned_int);
      subdomain_ids.assign(n_cells, 0);
      level_subdomain_ids.assign(n_cells, 0);

      refine_flags.assign(n_cells, 0u);
      refine_choice.assign(n_cells, 0u);
      coarsen_flags.assign(n_cells, false);

      parents.assign((n_cells + 1) / 2, -1);

      if (dim == spacedim - 1)
        direction_flags.assign(n_cells, true);

      cells.allocate(n_cells);

      face_orientations.reinit(orientation_needed ? n_cells : 0u,
                               faces_per_object);

      neighbors.assign(n_cells * faces_per_object, {-1, -1});

      reference_cell.assign(n_cells, ReferenceCells::Invalid<dim>);

      global_active_cell_indices.assign(n_cells, numbers::invalid_dof_index);
      global_level_cell_indices.assign(n_cells, numbers::invalid_dof_index);
    }



    template <int dim, int spacedim>
    void
    TriaLevel<dim, spacedim>::allocate_end(const std::size_t n_new_cells,
                                           const bool        orientation_needed,
                                           const bool        has_tetrahedra)
    {
      Assert(children_per_object == cells.children_per_object,
             ExcInternalError());
      Assert(faces_per_object == cells.faces_per_object, ExcInternalError());
      if (n_new_cells == 0)
        return;
      const auto total_cells = size() + n_new_cells;

      // note that all arrays should have equal sizes (checked by
      // @p{monitor_memory}
      refine_flags.reserve(total_cells);
      refine_flags.insert(refine_flags.end(),
                          total_cells - refine_flags.size(),
                          /*RefinementCase::no_refinement=*/0);

      if (has_tetrahedra)
        {
          refine_choice.reserve(total_cells);
          refine_choice.insert(
            refine_choice.end(),
            total_cells - refine_choice.size(),
            static_cast<std::uint8_t>(
              IsotropicRefinementChoice::isotropic_refinement));
        }

      coarsen_flags.reserve(total_cells);
      coarsen_flags.insert(coarsen_flags.end(),
                           total_cells - coarsen_flags.size(),
                           false);

      active_cell_indices.reserve(total_cells);
      active_cell_indices.insert(active_cell_indices.end(),
                                 total_cells - active_cell_indices.size(),
                                 numbers::invalid_unsigned_int);

      subdomain_ids.reserve(total_cells);
      subdomain_ids.insert(subdomain_ids.end(),
                           total_cells - subdomain_ids.size(),
                           0);

      level_subdomain_ids.reserve(total_cells);
      level_subdomain_ids.insert(level_subdomain_ids.end(),
                                 total_cells - level_subdomain_ids.size(),
                                 0);

      global_active_cell_indices.reserve(total_cells);
      global_active_cell_indices.insert(global_active_cell_indices.end(),
                                        total_cells -
                                          global_active_cell_indices.size(),
                                        numbers::invalid_dof_index);

      global_level_cell_indices.reserve(total_cells);
      global_level_cell_indices.insert(global_level_cell_indices.end(),
                                       total_cells -
                                         global_level_cell_indices.size(),
                                       numbers::invalid_dof_index);

      if (dim == spacedim - 1)
        {
          direction_flags.reserve(total_cells);
          direction_flags.insert(direction_flags.end(),
                                 total_cells - direction_flags.size(),
                                 true);
        }
      else
        direction_flags.clear();

      parents.reserve((total_cells + 1) / 2);
      parents.insert(parents.end(), (total_cells + 1) / 2 - parents.size(), -1);

      neighbors.reserve(total_cells * faces_per_object);
      neighbors.insert(neighbors.end(),
                       total_cells * faces_per_object - neighbors.size(),
                       std::make_pair(-1, -1));

      if (dim > 1)
        {
          face_orientations.resize(orientation_needed ? total_cells : 0);

          reference_cell.reserve(total_cells);
          reference_cell.insert(reference_cell.end(),
                                total_cells - reference_cell.size(),
                                ReferenceCells::get_hypercube<dim>());
        }
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
