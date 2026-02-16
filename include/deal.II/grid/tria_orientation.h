// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2023 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

#ifndef dealii_tria_orientation_h
#define dealii_tria_orientation_h

#include <deal.II/base/config.h>

#include <deal.II/base/types.h>
#include <deal.II/base/utilities.h>

#include <tuple>

DEAL_II_NAMESPACE_OPEN

namespace internal
{
  /**
   * Combine orientation flags.
   */
  inline types::geometric_orientation
  combined_face_orientation(const bool face_orientation,
                            const bool face_rotation,
                            const bool face_flip)
  {
    return (face_orientation ? 0 : 1) + (face_rotation ? 2 : 0) +
           (face_flip ? 4 : 0);
  }

  /**
   * Split up a combined orientation flag: orientation flag,
   * rotation flag, flip flag.
   */
  inline std::tuple<bool, bool, bool>
  split_face_orientation(
    const types::geometric_orientation combined_orientation)
  {
    Assert(combined_orientation != numbers::invalid_geometric_orientation,
           ExcMessage(
             "The provided orientation value is not valid. This typically "
             "means this value was initialized to an invalid value and not "
             "correctly set later."));
    return {!Utilities::get_bit(combined_orientation, 0),
            Utilities::get_bit(combined_orientation, 1),
            Utilities::get_bit(combined_orientation, 2)};
  }


} // namespace internal



DEAL_II_NAMESPACE_CLOSE

#endif
