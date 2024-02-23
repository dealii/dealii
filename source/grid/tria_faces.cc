// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2006 - 2022 by the deal.II authors
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
                MemoryConsumption::memory_consumption(quads_line_orientations) +
                MemoryConsumption::memory_consumption(quad_is_quadrilateral) +
                MemoryConsumption::memory_consumption(lines));

      return 0;
    }
  } // namespace TriangulationImplementation
} // namespace internal

DEAL_II_NAMESPACE_CLOSE
