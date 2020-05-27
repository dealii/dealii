// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2020 by the deal.II authors
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

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_objects.h>

#include <algorithm>
#include <functional>



DEAL_II_NAMESPACE_OPEN

namespace internal
{
  namespace TriangulationImplementation
  {
    void
    TriaObjects::reserve_space(const unsigned int new_objects_in_pairs,
                               const unsigned int new_objects_single)
    {
      if (structdim <= 2)
        {
          Assert(new_objects_in_pairs % 2 == 0, ExcInternalError());

          next_free_single               = 0;
          next_free_pair                 = 0;
          reverse_order_next_free_single = false;

          // count the number of objects, of unused single objects and of
          // unused pairs of objects
          unsigned int n_objects        = 0;
          unsigned int n_unused_pairs   = 0;
          unsigned int n_unused_singles = 0;
          for (unsigned int i = 0; i < used.size(); ++i)
            {
              if (used[i])
                ++n_objects;
              else if (i + 1 < used.size())
                {
                  if (used[i + 1])
                    {
                      ++n_unused_singles;
                      if (next_free_single == 0)
                        next_free_single = i;
                    }
                  else
                    {
                      ++n_unused_pairs;
                      if (next_free_pair == 0)
                        next_free_pair = i;
                      ++i;
                    }
                }
              else
                ++n_unused_singles;
            }
          Assert(n_objects + 2 * n_unused_pairs + n_unused_singles ==
                   used.size(),
                 ExcInternalError());

          // how many single objects are needed in addition to
          // n_unused_objects?
          const int additional_single_objects =
            new_objects_single - n_unused_singles;

          unsigned int new_size =
            used.size() + new_objects_in_pairs - 2 * n_unused_pairs;
          if (additional_single_objects > 0)
            new_size += additional_single_objects;

          // only allocate space if necessary
          if (new_size > this->n_objects())
            {
              unsigned int faces_per_cell        = 1;
              unsigned int max_children_per_cell = 1;

              if (this->structdim == 1)
                faces_per_cell = GeometryInfo<1>::faces_per_cell;
              else if (this->structdim == 2)
                faces_per_cell = GeometryInfo<2>::faces_per_cell;
              else if (this->structdim == 3)
                faces_per_cell = GeometryInfo<3>::faces_per_cell;
              else
                AssertThrow(false, ExcNotImplemented());

              if (this->structdim == 1)
                max_children_per_cell = GeometryInfo<1>::max_children_per_cell;
              else if (this->structdim == 2)
                max_children_per_cell = GeometryInfo<2>::max_children_per_cell;
              else if (this->structdim == 3)
                max_children_per_cell = GeometryInfo<3>::max_children_per_cell;
              else
                AssertThrow(false, ExcNotImplemented());

              cells.reserve(new_size * faces_per_cell);
              cells.insert(cells.end(),
                           (new_size - this->n_objects()) * faces_per_cell,
                           -1);

              used.reserve(new_size);
              used.insert(used.end(), new_size - used.size(), false);

              user_flags.reserve(new_size);
              user_flags.insert(user_flags.end(),
                                new_size - user_flags.size(),
                                false);

              const unsigned int factor = max_children_per_cell / 2;
              children.reserve(factor * new_size);
              children.insert(children.end(),
                              factor * new_size - children.size(),
                              -1);

              if (structdim > 1)
                {
                  refinement_cases.reserve(new_size);
                  refinement_cases.insert(refinement_cases.end(),
                                          new_size - refinement_cases.size(),
                                          /*RefinementCase::no_refinement=*/0);
                }

              // first reserve, then resize. Otherwise the std library can
              // decide to allocate more entries.
              boundary_or_material_id.reserve(new_size);
              boundary_or_material_id.resize(new_size);

              user_data.reserve(new_size);
              user_data.resize(new_size);

              manifold_id.reserve(new_size);
              manifold_id.insert(manifold_id.end(),
                                 new_size - manifold_id.size(),
                                 numbers::flat_manifold_id);
            }

          if (n_unused_singles == 0)
            {
              next_free_single               = new_size - 1;
              reverse_order_next_free_single = true;
            }
        }
      else
        {
          const unsigned int new_hexes = new_objects_in_pairs;

          const unsigned int new_size =
            new_hexes + std::count(used.begin(), used.end(), true);

          // see above...
          if (new_size > n_objects())
            {
              unsigned int faces_per_cell = 1;

              if (this->structdim == 1)
                faces_per_cell = GeometryInfo<1>::faces_per_cell;
              else if (this->structdim == 2)
                faces_per_cell = GeometryInfo<2>::faces_per_cell;
              else if (this->structdim == 3)
                faces_per_cell = GeometryInfo<3>::faces_per_cell;
              else
                AssertThrow(false, ExcNotImplemented());

              cells.reserve(new_size * faces_per_cell);
              cells.insert(cells.end(),
                           (new_size - n_objects()) * faces_per_cell,
                           -1);

              used.reserve(new_size);
              used.insert(used.end(), new_size - used.size(), false);

              user_flags.reserve(new_size);
              user_flags.insert(user_flags.end(),
                                new_size - user_flags.size(),
                                false);

              children.reserve(4 * new_size);
              children.insert(children.end(),
                              4 * new_size - children.size(),
                              -1);

              // for the following fields, we know exactly how many elements
              // we need, so first reserve then resize (resize itself, at least
              // with some compiler libraries, appears to round up the size it
              // actually reserves)
              boundary_or_material_id.reserve(new_size);
              boundary_or_material_id.resize(new_size);

              manifold_id.reserve(new_size);
              manifold_id.insert(manifold_id.end(),
                                 new_size - manifold_id.size(),
                                 numbers::flat_manifold_id);

              user_data.reserve(new_size);
              user_data.resize(new_size);

              refinement_cases.reserve(new_size);
              refinement_cases.insert(refinement_cases.end(),
                                      new_size - refinement_cases.size(),
                                      /*RefinementCase::no_refinement=*/0);
            }
          next_free_single = next_free_pair = 0;
        }
    }


    template <int dim, int spacedim>
    typename dealii::Triangulation<dim, spacedim>::raw_hex_iterator
    TriaObjects::next_free_hex(const dealii::Triangulation<dim, spacedim> &tria,
                               const unsigned int level)
    {
      AssertDimension(this->structdim, 3);

      // TODO: Think of a way to ensure that we are using the correct
      // triangulation, i.e. the one containing *this.

      int pos = next_free_pair, last = used.size() - 1;
      for (; pos < last; ++pos)
        if (!used[pos])
          {
            // this should be a pair slot
            Assert(!used[pos + 1], ExcInternalError());
            break;
          }
      if (pos >= last)
        // no free slot
        return tria.end_hex();
      else
        next_free_pair = pos + 2;

      return
        typename dealii::Triangulation<dim, spacedim>::raw_hex_iterator(&tria,
                                                                        level,
                                                                        pos);
    }


    void
    TriaObjects::monitor_memory(const unsigned int) const
    {
      if (this->structdim == 1)
        {
          Assert(n_objects() == used.size(),
                 ExcMemoryInexact(n_objects(), used.size()));
          Assert(n_objects() == user_flags.size(),
                 ExcMemoryInexact(n_objects(), user_flags.size()));
          Assert(n_objects() == children.size(),
                 ExcMemoryInexact(n_objects(), children.size()));
          Assert(n_objects() == boundary_or_material_id.size(),
                 ExcMemoryInexact(n_objects(), boundary_or_material_id.size()));
          Assert(n_objects() == manifold_id.size(),
                 ExcMemoryInexact(n_objects(), manifold_id.size()));
          Assert(n_objects() == user_data.size(),
                 ExcMemoryInexact(n_objects(), user_data.size()));
        }
      else if (this->structdim == 2)
        {
          Assert(n_objects() == used.size(),
                 ExcMemoryInexact(n_objects(), used.size()));
          Assert(n_objects() == user_flags.size(),
                 ExcMemoryInexact(n_objects(), user_flags.size()));
          Assert(2 * n_objects() == children.size(),
                 ExcMemoryInexact(n_objects(), children.size()));
          Assert(n_objects() == refinement_cases.size(),
                 ExcMemoryInexact(n_objects(), refinement_cases.size()));
          Assert(n_objects() == boundary_or_material_id.size(),
                 ExcMemoryInexact(n_objects(), boundary_or_material_id.size()));
          Assert(n_objects() == manifold_id.size(),
                 ExcMemoryInexact(n_objects(), manifold_id.size()));
          Assert(n_objects() == user_data.size(),
                 ExcMemoryInexact(n_objects(), user_data.size()));
        }
      else if (this->structdim == 3)
        {
          Assert(n_objects() == used.size(),
                 ExcMemoryInexact(n_objects(), used.size()));
          Assert(n_objects() == user_flags.size(),
                 ExcMemoryInexact(n_objects(), user_flags.size()));
          Assert(4 * n_objects() == children.size(),
                 ExcMemoryInexact(n_objects(), children.size()));
          Assert(n_objects() == boundary_or_material_id.size(),
                 ExcMemoryInexact(n_objects(), boundary_or_material_id.size()));
          Assert(n_objects() == manifold_id.size(),
                 ExcMemoryInexact(n_objects(), manifold_id.size()));
          Assert(n_objects() == user_data.size(),
                 ExcMemoryInexact(n_objects(), user_data.size()));
        }
    }


    std::size_t
    TriaObjects::memory_consumption() const
    {
      return (MemoryConsumption::memory_consumption(cells) +
              MemoryConsumption::memory_consumption(children) +
              MemoryConsumption::memory_consumption(used) +
              MemoryConsumption::memory_consumption(user_flags) +
              MemoryConsumption::memory_consumption(boundary_or_material_id) +
              MemoryConsumption::memory_consumption(manifold_id) +
              MemoryConsumption::memory_consumption(refinement_cases) +
              user_data.capacity() * sizeof(UserData) + sizeof(user_data));
    }



    // explicit instantiations
#ifndef DOXYGEN
#  include "tria_objects.inst"
#endif
  } // namespace TriangulationImplementation
} // namespace internal

DEAL_II_NAMESPACE_CLOSE
