// ---------------------------------------------------------------------
//
// Copyright (C) 2020 by the deal.II authors
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

#ifndef dealii_tria_referce_cell_h
#define dealii_tria_referce_cell_h

#include <deal.II/base/config.h>


DEAL_II_NAMESPACE_OPEN


/**
 * A namespace for reference cells.
 */
namespace ReferenceCell
{
  /**
   * Supported reference cell types.
   */
  enum class Type : std::uint8_t
  {
    Vertex  = 0,
    Line    = 1,
    Tri     = 2,
    Quad    = 3,
    Tet     = 4,
    Pyramid = 5,
    Wedge   = 6,
    Hex     = 7,
    Invalid = static_cast<std::uint8_t>(-1)
  };
} // namespace ReferenceCell


DEAL_II_NAMESPACE_CLOSE

#endif
