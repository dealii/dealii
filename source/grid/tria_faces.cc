// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2014 by the deal.II authors
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
#include <deal.II/grid/tria_faces.h>

DEAL_II_NAMESPACE_OPEN


namespace internal
{
  namespace Triangulation
  {

    std::size_t
    TriaFaces<1>::memory_consumption () const
    {
      return 0;
    }


    std::size_t
    TriaFaces<2>::memory_consumption () const
    {
      return MemoryConsumption::memory_consumption (lines);
    }


    std::size_t
    TriaFaces<3>::memory_consumption () const
    {
      return (MemoryConsumption::memory_consumption (quads) +
              MemoryConsumption::memory_consumption (lines) );
    }
  }
}

DEAL_II_NAMESPACE_CLOSE

