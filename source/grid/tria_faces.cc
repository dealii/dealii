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
  } // namespace TriangulationImplementation
} // namespace internal

DEAL_II_NAMESPACE_CLOSE
