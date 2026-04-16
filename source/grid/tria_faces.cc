// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2006 - 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


#include <deal.II/base/memory_consumption.h>

#include <deal.II/grid/tria_faces.h>

DEAL_II_NAMESPACE_OPEN


namespace internal
{
  namespace TriangulationImplementation
  {
    template <int dim>
    TriaFaces<dim>::TriaFaces(const unsigned int max_children_per_quad,
                              const unsigned int max_lines_per_quad)
      : max_lines_per_quad(max_lines_per_quad)
      , quads(2, max_children_per_quad, max_lines_per_quad)
      , lines(1,
              ReferenceCells::max_n_children<1>(),
              ReferenceCells::max_n_faces<1>())
    {}



    template <int dim>
    void
    TriaFaces<dim>::allocate(const std::size_t n_lines,
                             const std::size_t n_quads)
    {
      if (dim > 1)
        {
          lines.allocate(n_lines);
        }

      if (dim == 3)
        {
          quads.allocate(n_quads);
          quad_is_quadrilateral.assign(n_quads, true);
          quads_line_orientations.assign(n_quads * max_lines_per_quad, true);
        }
    }

    template <int dim>
    void
    TriaFaces<dim>::allocate_end(const std::size_t new_quads_in_pairs,
                                 const std::size_t new_quads_single)
    {
      if (dim < 3)
        return;

      Assert(new_quads_in_pairs % 2 == 0, ExcInternalError());

      std::size_t next_free_single = 0;
      std::size_t next_free_pair   = 0;

      // count the number of objects, of unused single objects and of
      // unused pairs of objects
      std::size_t    n_quads          = 0;
      std::size_t    n_unused_pairs   = 0;
      std::ptrdiff_t n_unused_singles = 0;
      for (std::size_t i = 0; i < quads.used.size(); ++i)
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
      const auto additional_single_quads =
        std::ptrdiff_t(new_quads_single) - n_unused_singles;

      Assert(quads.used.size() + new_quads_in_pairs >= 2 * n_unused_pairs,
             ExcInternalError());
      auto new_size =
        quads.used.size() + new_quads_in_pairs - 2 * n_unused_pairs;
      if (additional_single_quads > 0)
        new_size += additional_single_quads;

      // see above...
      if (new_size > quads.n_objects())
        {
          // reserve the field of the derived class
          quads_line_orientations.resize(new_size * quads.faces_per_object,
                                         true);

          quad_is_quadrilateral.reserve(new_size);
          quad_is_quadrilateral.insert(quad_is_quadrilateral.end(),
                                       new_size - quad_is_quadrilateral.size(),
                                       true);
        }
    }



    template <int dim>
    std::size_t
    TriaFaces<dim>::memory_consumption() const
    {
      if (dim == 2)
        return MemoryConsumption::memory_consumption(lines);
      if (dim == 3)
        return (MemoryConsumption::memory_consumption(quads) +
                MemoryConsumption::memory_consumption(quads_line_orientations) +
                MemoryConsumption::memory_consumption(quad_is_quadrilateral) +
                MemoryConsumption::memory_consumption(lines));

      return 0;
    }


    template class TriaFaces<1>;
    template class TriaFaces<2>;
    template class TriaFaces<3>;
  } // namespace TriangulationImplementation
} // namespace internal

DEAL_II_NAMESPACE_CLOSE
