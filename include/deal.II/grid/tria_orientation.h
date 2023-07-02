// ---------------------------------------------------------------------
//
// Copyright (C) 2023 by the deal.II authors
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

#ifndef dealii_tria_orientation_h
#define dealii_tria_orientation_h

#include <deal.II/base/config.h>

#include <deal.II/base/utilities.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN

namespace internal
{
  /**
   * Combine orientation flags.
   */
  inline unsigned char
  combined_face_orientation(const bool face_orientation,
                            const bool face_rotation,
                            const bool face_flip)
  {
    return (face_orientation ? 1 : 0) + (face_rotation ? 2 : 0) +
           (face_flip ? 4 : 0);
  }

  /**
   * Split up a combined orientation flag: orientation flag,
   * rotation flag, flip flag.
   */
  inline std::tuple<bool, bool, bool>
  split_face_orientation(const unsigned char combined_face_orientation)
  {
    return {Utilities::get_bit(combined_face_orientation, 0),
            Utilities::get_bit(combined_face_orientation, 1),
            Utilities::get_bit(combined_face_orientation, 2)};
  }


} // namespace internal



DEAL_II_NAMESPACE_CLOSE

#endif
