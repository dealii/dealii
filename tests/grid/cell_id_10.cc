// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


// Verify that CellId works with very large coarse cell ids in 64-bit mode.
// Previously, we only used an unsigned int to encode the coarse cell index,
// which isn't large enough for all possible coarse cell indices.

#include <deal.II/grid/cell_id.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <sstream>

#include "../tests.h"

void
test(const std::uint8_t n_refined_levels)
{
  auto test_id = [&](const CellId cell_id_1) {
    deallog << cell_id_1 << std::endl;
    const CellId cell_id_2(cell_id_1.to_binary<3>());
    deallog << cell_id_2 << std::endl;
    AssertThrow(cell_id_1 == cell_id_2, ExcInternalError());

    std::ostringstream out_stream;
    out_stream << cell_id_1;
    std::istringstream in_stream(out_stream.str());
    CellId             cell_id_3;
    in_stream >> cell_id_3;
    AssertThrow(cell_id_1 == cell_id_3, ExcInternalError());
  };

  // 2d:
  {
    const auto coarse_index =
      Utilities::fixed_power<60, types::coarse_cell_id>(2) +
      Utilities::fixed_power<10, types::coarse_cell_id>(2) + 1;

    std::vector<std::uint8_t> child_indices;
    for (std::uint8_t i = 0; i < n_refined_levels; ++i)
      child_indices.push_back((3 * i + 1) % 4);

    const CellId cell_id_1(coarse_index, child_indices);
    test_id(cell_id_1);
  }

  // 3d (anything but Pyramid):
  {
    const auto coarse_index =
      Utilities::fixed_power<50, types::coarse_cell_id>(2) +
      Utilities::fixed_power<10, types::coarse_cell_id>(2) + 1;

    std::vector<std::uint8_t> child_indices;
    for (std::uint8_t i = 0; i < n_refined_levels; ++i)
      child_indices.push_back((3 * i + 1) % 8);

    const CellId cell_id_1(coarse_index, child_indices);
    test_id(cell_id_1);
  }

  // 3d, based on Pyramid:
#if 0
  {
    const auto coarse_index =
      Utilities::fixed_power<50, types::coarse_cell_id>(2) +
      Utilities::fixed_power<10, types::coarse_cell_id>(2) + 1;

    std::vector<std::uint8_t> child_indices;
    for (std::uint8_t i = 0; i < n_refined_levels; ++i)
      child_indices.push_back((3 * i + 1) % 10);

    const CellId cell_id_1(coarse_index, child_indices);
    test_id(cell_id_1);
  }
#endif
}


int
main()
{
  initlog();

  deallog << "one level:" << std::endl;
  test(0);
  deallog << "two levels:" << std::endl;
  test(1);
  deallog << "20 levels:" << std::endl;
  test(19);
  deallog << "30 levels:" << std::endl;
  test(29);
}
