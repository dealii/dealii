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
    TriaFaces::TriaFaces(const unsigned int dim,
                         const unsigned int max_children_per_quad,
                         const unsigned int max_lines_per_quad)
      : dim(dim)
      , quads(2, max_children_per_quad, max_lines_per_quad)
      , lines(1,
              ReferenceCells::max_n_children<1>(),
              ReferenceCells::max_n_faces<1>())
    {}

    std::size_t
    TriaFaces::memory_consumption() const
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
  } // namespace TriangulationImplementation
} // namespace internal

DEAL_II_NAMESPACE_CLOSE
