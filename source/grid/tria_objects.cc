// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2006 - 2025 by the deal.II authors
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
#  include "grid/tria_objects.inst"
#endif
  } // namespace TriangulationImplementation
} // namespace internal

DEAL_II_NAMESPACE_CLOSE
