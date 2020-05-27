// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2018 by the deal.II authors
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

#include <deal.II/grid/tria_faces.h>

DEAL_II_NAMESPACE_OPEN


namespace internal
{
  namespace TriangulationImplementation
  {
    TriaFaces::TriaFaces(const unsigned int dim)
      : dim(dim)
      , quads(2)
      , lines(1)
    {}

    std::size_t
    TriaFaces::memory_consumption() const
    {
      if (dim == 2)
        return MemoryConsumption::memory_consumption(lines);
      if (dim == 3)
        return (MemoryConsumption::memory_consumption(quads) +
                MemoryConsumption::memory_consumption(lines));

      return 0;
    }

    void
    TriaFaces::reserve_space(const unsigned int new_quads_in_pairs,
                             const unsigned int new_quads_single)
    {
      AssertDimension(this->dim, 3);

      Assert(new_quads_in_pairs % 2 == 0, ExcInternalError());

      unsigned int next_free_single = 0;
      unsigned int next_free_pair   = 0;

      // count the number of objects, of unused single objects and of
      // unused pairs of objects
      unsigned int n_quads          = 0;
      unsigned int n_unused_pairs   = 0;
      unsigned int n_unused_singles = 0;
      for (unsigned int i = 0; i < quads.used.size(); ++i)
        {
          if (quads.used[i])
            ++n_quads;
          else if (i + 1 < quads.used.size())
            {
              if (quads.used[i + 1])
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
      Assert(n_quads + 2 * n_unused_pairs + n_unused_singles ==
               quads.used.size(),
             ExcInternalError());

      // how many single quads are needed in addition to n_unused_quads?
      const int additional_single_quads = new_quads_single - n_unused_singles;

      unsigned int new_size =
        quads.used.size() + new_quads_in_pairs - 2 * n_unused_pairs;
      if (additional_single_quads > 0)
        new_size += additional_single_quads;

      // see above...
      if (new_size > quads.n_objects())
        {
          // reserve the field of the derived class
          quads_line_orientations.reserve(new_size *
                                          GeometryInfo<2>::lines_per_cell);
          quads_line_orientations.insert(quads_line_orientations.end(),
                                         new_size *
                                             GeometryInfo<2>::lines_per_cell -
                                           quads_line_orientations.size(),
                                         true);
        }
    }
  } // namespace TriangulationImplementation
} // namespace internal

DEAL_II_NAMESPACE_CLOSE
