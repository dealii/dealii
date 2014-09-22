// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2013 by the deal.II authors
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
#include <deal.II/grid/tria_levels.h>

DEAL_II_NAMESPACE_OPEN


namespace internal
{
  namespace Triangulation
  {
    template <int dim>
    void
    TriaLevel<dim>::reserve_space (const unsigned int total_cells,
                                   const unsigned int dimension,
                                   const unsigned int space_dimension)
    {
      // we need space for total_cells cells. Maybe we have more already
      // with those cells which are unused, so only allocate new space if
      // needed.
      //
      // note that all arrays should have equal sizes (checked by
      // @p{monitor_memory}
      if (total_cells > refine_flags.size())
        {
          refine_flags.reserve (total_cells);
          refine_flags.insert (refine_flags.end(),
                               total_cells - refine_flags.size(),
                               RefinementCase<dim>::no_refinement);

          coarsen_flags.reserve (total_cells);
          coarsen_flags.insert (coarsen_flags.end(),
                                total_cells - coarsen_flags.size(),
                                false);

          subdomain_ids.reserve (total_cells);
          subdomain_ids.insert (subdomain_ids.end(),
                                total_cells - subdomain_ids.size(),
                                0);

          level_subdomain_ids.reserve (total_cells);
          level_subdomain_ids.insert (level_subdomain_ids.end(),
                                      total_cells - level_subdomain_ids.size(),
                                      0);

          if (dimension < space_dimension)
            {
              direction_flags.reserve (total_cells);
              direction_flags.insert (direction_flags.end(),
                                      total_cells - direction_flags.size(),
                                      true);
            }
          else
            direction_flags.clear ();

          parents.reserve ((int) (total_cells + 1) / 2);
          parents.insert (parents.end (),
                          (total_cells + 1) / 2 - parents.size (),
                          -1);

          neighbors.reserve (total_cells*(2*dimension));
          neighbors.insert (neighbors.end(),
                            total_cells*(2*dimension) - neighbors.size(),
                            std::make_pair(-1,-1));
        };
    }


    template <int dim>
    void
    TriaLevel<dim>::monitor_memory (const unsigned int true_dimension) const
    {
      Assert (2*true_dimension*refine_flags.size() == neighbors.size(),
              ExcMemoryInexact (refine_flags.size(), neighbors.size()));
      Assert (2*true_dimension*coarsen_flags.size() == neighbors.size(),
              ExcMemoryInexact (coarsen_flags.size(), neighbors.size()));
    }


    template <int dim>
    std::size_t
    TriaLevel<dim>::memory_consumption () const
    {
      return (MemoryConsumption::memory_consumption (refine_flags) +
              MemoryConsumption::memory_consumption (coarsen_flags) +
              MemoryConsumption::memory_consumption (neighbors) +
              MemoryConsumption::memory_consumption (subdomain_ids) +
              MemoryConsumption::memory_consumption (level_subdomain_ids) +
              MemoryConsumption::memory_consumption (parents) +
              MemoryConsumption::memory_consumption (direction_flags) +
              MemoryConsumption::memory_consumption (cells));
    }

// This specialization should be only temporary, until the TriaObjects
// classes are straightened out.

    void
    TriaLevel<3>::reserve_space (const unsigned int total_cells,
                                 const unsigned int dimension,
                                 const unsigned int space_dimension)
    {
      // we need space for total_cells
      // cells. Maybe we have more already
      // with those cells which are unused,
      // so only allocate new space if needed.
      //
      // note that all arrays should have equal
      // sizes (checked by @p{monitor_memory}
      if (total_cells > refine_flags.size())
        {
          refine_flags.reserve (total_cells);
          refine_flags.insert (refine_flags.end(),
                               total_cells - refine_flags.size(),
                               RefinementCase<3>::no_refinement);

          coarsen_flags.reserve (total_cells);
          coarsen_flags.insert (coarsen_flags.end(),
                                total_cells - coarsen_flags.size(),
                                false);

          subdomain_ids.reserve (total_cells);
          subdomain_ids.insert (subdomain_ids.end(),
                                total_cells - subdomain_ids.size(),
                                0);

          level_subdomain_ids.reserve (total_cells);
          level_subdomain_ids.insert (level_subdomain_ids.end(),
                                      total_cells - level_subdomain_ids.size(),
                                      0);

          if (dimension < space_dimension)
            {
              direction_flags.reserve (total_cells);
              direction_flags.insert (direction_flags.end(),
                                      total_cells - direction_flags.size(),
                                      true);
            }
          else
            direction_flags.clear ();

          parents.reserve ((int) (total_cells + 1) / 2);
          parents.insert (parents.end (),
                          (total_cells + 1) / 2 - parents.size (),
                          -1);

          neighbors.reserve (total_cells*(2*dimension));
          neighbors.insert (neighbors.end(),
                            total_cells*(2*dimension) - neighbors.size(),
                            std::make_pair(-1,-1));
        };
    }


    void
    TriaLevel<3>::monitor_memory (const unsigned int true_dimension) const
    {
      Assert (2*true_dimension*refine_flags.size() == neighbors.size(),
              ExcMemoryInexact (refine_flags.size(), neighbors.size()));
      Assert (2*true_dimension*coarsen_flags.size() == neighbors.size(),
              ExcMemoryInexact (coarsen_flags.size(), neighbors.size()));
    }


    std::size_t
    TriaLevel<3>::memory_consumption () const
    {
      return (MemoryConsumption::memory_consumption (refine_flags) +
              MemoryConsumption::memory_consumption (coarsen_flags) +
              MemoryConsumption::memory_consumption (neighbors) +
              MemoryConsumption::memory_consumption (subdomain_ids) +
              MemoryConsumption::memory_consumption (parents) +
              MemoryConsumption::memory_consumption (direction_flags) +
              MemoryConsumption::memory_consumption (cells));
    }
  }
}


template class internal::Triangulation::TriaLevel<1>;
template class internal::Triangulation::TriaLevel<2>;

DEAL_II_NAMESPACE_CLOSE

