// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2024 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test that we get the correct reversed orientations (i.e., the composition of
// permutations is the identity operation)

#include <deal.II/grid/reference_cell.h>
#include <deal.II/grid/tria_orientation.h>

#include "../tests.h"

int
main()
{
  initlog();

  deallog << "lines" << std::endl;
  {
    const auto reference_cell = ReferenceCells::Line;
    for (types::geometric_orientation o = 0; o < 2; ++o)
      {
        deallog << "o = " << int(o) << std::endl;
        const std::array<unsigned int, 2> vs{{0u, 1u}};

        const auto result1 = reference_cell.permute_by_combined_orientation(
          make_const_array_view(vs), o);

        ArrayView<const unsigned int> view1(result1.data(), result1.size());
        const auto result2 = reference_cell.permute_by_combined_orientation(
          view1, reference_cell.get_inverse_combined_orientation(o));

        for (const auto &v : result2)
          deallog << "  " << v << std::endl;
      }
  }

  deallog << "triangles" << std::endl;
  {
    const auto reference_cell = ReferenceCells::Triangle;
    for (types::geometric_orientation o = 0; o < 6; ++o)
      {
        deallog << "o = " << int(o) << std::endl;
        const std::array<unsigned int, 3> vs{{0u, 1u, 2u}};

        const auto result1 = reference_cell.permute_by_combined_orientation(
          make_const_array_view(vs), o);

        ArrayView<const unsigned int> view1(result1.data(), result1.size());
        const auto result2 = reference_cell.permute_by_combined_orientation(
          view1, reference_cell.get_inverse_combined_orientation(o));

        for (const auto &v : result2)
          deallog << "  " << v << std::endl;
      }
  }

  deallog << "quadrilaterals" << std::endl;
  {
    const auto reference_cell = ReferenceCells::Quadrilateral;
    for (types::geometric_orientation o = 0; o < 8; ++o)
      {
        deallog << "o = " << int(o) << std::endl;
        const std::array<unsigned int, 4> vs{{0u, 1u, 2u, 3u}};

        const auto result1 = reference_cell.permute_by_combined_orientation(
          make_const_array_view(vs), o);

        ArrayView<const unsigned int> view1(result1.data(), result1.size());
        const auto result2 = reference_cell.permute_by_combined_orientation(
          view1, reference_cell.get_inverse_combined_orientation(o));

        for (const auto &v : result2)
          deallog << "  " << v << std::endl;
      }
  }

  // Verify that the manual version created the same results.
  deallog << "quadrilaterals (manual)" << std::endl;
  {
    const auto reference_cell = ReferenceCells::Quadrilateral;
    for (types::geometric_orientation o = 0; o < 8; ++o)
      {
        deallog << "o = " << int(o) << std::endl;
        const std::array<unsigned int, 4> vs{{0u, 1u, 2u, 3u}};

        const auto result1 = reference_cell.permute_by_combined_orientation(
          make_const_array_view(vs), o);

        ArrayView<const unsigned int> view1(result1.data(), result1.size());

        auto [orientation, rotation, flip] =
          internal::split_face_orientation(o);
        flip               = orientation ? rotation ^ flip : flip;
        const auto result2 = reference_cell.permute_by_combined_orientation(
          view1,
          internal::combined_face_orientation(orientation, rotation, flip));

        for (const auto &v : result2)
          deallog << "  " << v << std::endl;
      }
  }
}
